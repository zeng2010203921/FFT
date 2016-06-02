#define _GNU_SOURCE 
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <pthread.h>
#include <error.h>
#include <unistd.h>

#define BM 0x4D42					//位图bmp文件的ASCII码
#define PI 3.1415926535898
#define Unit16 unsigned int 

int cpu_nub=0;						//获得系统中cpu的个数

unsigned int OffSet=0;				//位图文件头、信息头的偏移量
long bmpWidth=0;					//图像宽度
long bmpHeight=0;					//图像高度
unsigned int bmpBitCount=0;			//图像的位深度
unsigned char  *temp_bmp_arry;		//存放像素点的一维数组

long int line_data;					//每行应读入的数据的长度
long int Sample_Num;				//FFT变换数据点的个数
long int array_line_data;			//图像长*宽的长度
long int pixel_length;				//图像所包含的总的像素点的数量

void init_temp_array(double *p)     //对辅助一维数组进行初始化
{
	line_data=(bmpWidth*bmpBitCount+31)/8;
	line_data=line_data/4*4;
	temp_bmp_arry=(unsigned char *)malloc(sizeof(unsigned char)*line_data);//用于存取图像每一行的数据
}

/*快速傅里叶变换相关函数*/

struct FFT_DATA						//复数结构体
{
	float data_R;					//实部
	float data_I;					//虚部
};

struct WN
{
	float sin;
	float cos;
};

struct FFT_DATA *fft_r_data;
struct FFT_DATA *fft_g_data;
struct FFT_DATA *fft_b_data;
struct WN *Wn;

void my_init()
{
	fft_r_data=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*array_line_data);
	fft_g_data=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*array_line_data);
	fft_b_data=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*array_line_data);
	Wn=(struct WN *)malloc(sizeof(struct WN)*Sample_Num);
}

void my_free()
{
	free(fft_r_data);
	free(fft_g_data);
	free(fft_b_data);
	free(Wn);
}

void InitForFFT()						
{
	int i;
	for(i=0;i<Sample_Num;i++)			//(2*PI*n/N)——与公式相比，缺少一个k
	{
		Wn[i].sin=sin(PI*2*i/Sample_Num);
		Wn[i].cos=cos(PI*2*i/Sample_Num);
	}
}

void FFT_N(struct FFT_DATA *fft_data) //R G B分别进行FFT变换
{
	Unit16 xx,Num;
	Unit16 i,j,k,b,p,L;
	float TR,TI,temp;
	Num=(Unit16)(log(Sample_Num)/log(2)+0.0001);	//num=log(Sample_Num)

	//数据的倒序与原位
	for(i=0;i<Sample_Num;i++)
	{
		xx=0;
		for(j=0;j<Num;j++)
		{
			if(i&(1<<j))
			{
				xx|=(1<<(Num-1-j));
			}
		}

		fft_data[xx].data_I=fft_data[i].data_R;
	}

	for(i=0;i<Sample_Num;i++)
	{
		fft_data[i].data_R=fft_data[i].data_I;
		fft_data[i].data_I=0.0;
	}
	/*
	   蝶形运算:
	   L是第L级蝶形运算，与Num的值有关，总共是Num级蝶形运算，Num=log(N);
	   b是第L级做蝶形运算的两个点之间的距离，p是蝶形运算的旋转因子Wnp，根据旋转因子的对称性可知，旋转		　因子只有N/2个值;
	   由于第K个点与第K+b个点进行运算的结果仍然要保留在第K个点的位置，所以使用中间变量TR,TI,temp，防		  止结果被覆盖;

	   FFT算法旋转因子规律:
	   第一级蝶形系数p为0,蝶形节点的距离为1;
	   第二级蝶形系数p为0.N/4，蝶形节点的距离为2;
	   第三级蝶形系数p为0,N/8，2N/8,3N/8，蝶形节点的距离为4;
	   ''''''
	   */

	for(L=1;L<=Num;L++)//L表示第几级
	{
		b=pow(2,L-1);//计算第L级蝶形运算的距离
		for(j=0;j<=b-1;j++)//根据b的值，求蝶形系数p的范围
		{
			p=pow(2,Num-L)*j;			  //
			for(k=j;k<Sample_Num;k=k+2*b)//计算每一级蝶形运算的旋转因子,第k与第k+b距离的两个点进行蝶形											运算
			{
				TR=fft_data[k].data_R;
				TI=fft_data[k].data_I;
				temp=fft_data[k+b].data_R;

				fft_data[k].data_R=fft_data[k].data_R+fft_data[k+b].data_R*Wn[p].cos+fft_data[k+b].data_I*Wn[p].sin;

				fft_data[k].data_I=fft_data[k].data_I-fft_data[k+b].data_R*Wn[p].sin+fft_data[k+b].data_I*Wn[p].cos;

				fft_data[k+b].data_R=TR-fft_data[k+b].data_R*Wn[p].cos-fft_data[k+b].data_I*Wn[p].sin;

				fft_data[k+b].data_I=TI+temp*Wn[p].sin-fft_data[k+b].data_I*Wn[p].cos;
			}
		}
	}
}

struct mypara
{
	int thread_id;	//线程号
	unsigned char *data;
};

void *thread_fun_r(void *arg)
{

	struct mypara *pstru;
	pstru=(struct mypara *)arg;
	int th_id=(*pstru).thread_id;
	unsigned char *fp_r_data=(*pstru).data;

	cpu_set_t mask;		//创建cpu核的集合
	cpu_set_t get;		//获取在集合中的cpu
	CPU_ZERO(&mask);	//置空
	CPU_SET(th_id,&mask);
	if(sched_setaffinity(0,sizeof(mask),&mask)==-1)
	{
		printf("警告:不能设置线程的亲和性,继续...\n");
	}
	while(1)
	{
		CPU_ZERO(&get);
		if(sched_getaffinity(0,sizeof(get),&get)==-1)
		{
			printf("警告:不能获得线程的亲和性,继续...\n");
		}
		int t;
		for(t=0;t<cpu_nub;t++)
		{
			if(CPU_ISSET(t,&get))
			{
		 		printf("线程%d在第%dCPU上运行\n",th_id,t);
				break;
			}
		}
		break;
	}
	printf("R设置完毕\n");

	int u,v;
	struct FFT_DATA *r_arry=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*bmpHeight);
	for(u=0;u<bmpWidth;u++)
	{
		for(v=0;v<bmpHeight;v++)
		{
			r_arry[v].data_R=fp_r_data[v*bmpWidth+u];
			r_arry[v].data_I=0.0;
		}
		FFT_N(r_arry);
		for(v=0;v<bmpHeight;v++)
		{
			fft_r_data[v*bmpWidth+u].data_R=r_arry[v].data_R;
			fft_r_data[v*bmpWidth+u].data_I=r_arry[v].data_I;
		}
	}
	for(v=0;v<bmpHeight;v++)
	{
		FFT_N(fft_r_data+v*bmpWidth);
	}

	//数据变换后存在于R\G\B各自的一维输出数组中
	int r1;
	for(int j=0;j<bmpHeight;j++)
	{
		for(int i=0;i<bmpWidth;i++)
		{
			double temp1=fft_r_data[i+j*bmpWidth].data_R*fft_r_data[i+j*bmpWidth].data_R+fft_r_data[i+j*bmpWidth].data_I*fft_r_data[i+j*bmpWidth].data_I;
			r1=(int)(sqrt(temp1)/100);
			int u1,v1;
			if(i<bmpWidth/2)
				u1=i+bmpWidth/2;
			else
				u1=i-bmpWidth/2;
			if(j<bmpHeight/2)
				v1=j+bmpHeight/2;
			else
				v1=j-bmpHeight/2;
			if(r1>255)
				r1=255;
			if(r1<0)
				r1=0;
			fp_r_data[u1+v1*bmpWidth]=r1;
		}
	}
	//pthread_exit(1);
	pthread_exit(NULL);
}

void *thread_fun_g(void *arg)
{

	struct mypara *pstru;
	pstru=(struct mypara *)arg;
	int th_id=(*pstru).thread_id;
	unsigned char *fp_g_data=(*pstru).data;

	cpu_set_t mask;		//创建cpu核的集合
	cpu_set_t get;		//获取在集合中的cpu
	CPU_ZERO(&mask);	//置空
	CPU_SET(th_id,&mask);
	if(sched_setaffinity(0,sizeof(mask),&mask)==-1)
	{
		printf("警告:不能设置线程的亲和性,继续...\n");
	}
	while(1)
	{
		CPU_ZERO(&get);
		if(sched_getaffinity(0,sizeof(get),&get)==-1)
		{
			printf("警告:不能获得线程的亲和性,继续...\n");
		}
		int t;
		for(t=0;t<cpu_nub;t++)
		{
			if(CPU_ISSET(t,&get))
			{
				printf("线程%d在第%dCPU上运行\n",th_id,t);
				break;
			}
		}
		break;
	}
	printf("G设置完毕\n");

	struct FFT_DATA *g_arry=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*bmpHeight);
	int u,v;
	for(u=0;u<bmpWidth;u++)
	{
		for(v=0;v<bmpHeight;v++)
		{
			g_arry[v].data_R=fp_g_data[v*bmpWidth+u];
			g_arry[v].data_I=0.0;
		}
		FFT_N(g_arry);
		for(v=0;v<bmpHeight;v++)
		{
			fft_g_data[v*bmpWidth+u].data_R=g_arry[v].data_R;
			fft_g_data[v*bmpWidth+u].data_I=g_arry[v].data_I;
		}
	}
	for(v=0;v<bmpHeight;v++)
	{
		FFT_N(fft_g_data+v*bmpWidth);
	}

	//数据变换后存在于R\G\B各自的一维输出数组中
	int r2;
	for(int j=0;j<bmpHeight;j++)
	{
		for(int i=0;i<bmpWidth;i++)
		{
			double temp2=fft_g_data[i+j*bmpWidth].data_R*fft_g_data[i+j*bmpWidth].data_R+fft_g_data[i+j*bmpWidth].data_I*fft_g_data[i+j*bmpWidth].data_I;
			r2=(int)(sqrt(temp2)/100);
			int u1,v1;
			if(i<bmpWidth/2)
				u1=i+bmpWidth/2;
			else
				u1=i-bmpWidth/2;
			if(j<bmpHeight/2)
				v1=j+bmpHeight/2;
			else
				v1=j-bmpHeight/2;
			if(r2>255)
				r2=255;
			if(r2<0)			
				r2=0;
			fp_g_data[u1+v1*bmpWidth]=r2;
		}
	}
	//pthread_exit(2);
	pthread_exit(NULL);
}

void *thread_fun_b(void *arg)
{
	struct mypara *pstru;
	pstru=(struct mypara *)arg;
	int th_id=(*pstru).thread_id;
	unsigned char *fp_b_data=(*pstru).data;

	cpu_set_t mask;		//创建cpu核的集合
	cpu_set_t get;		//获取在集合中的cpu
	CPU_ZERO(&mask);	//置空
	CPU_SET(th_id,&mask);
	if(sched_setaffinity(0,sizeof(mask),&mask)==-1)
	{
		printf("警告:不能设置线程的亲和性,继续...\n");
	}
	while(1)
	{
		CPU_ZERO(&get);
		if(sched_getaffinity(0,sizeof(get),&get)==-1)
		{
			printf("警告:不能获得线程的亲和性,继续...\n");
		}
		int t;
		for(t=0;t<cpu_nub;t++)
		{
			if(CPU_ISSET(t,&get))
			{
				printf("线程%d在第%dCPU上运行\n",th_id,t);
				break;
			}
		}
		break;
	}
	printf("B设置完毕\n");

	struct FFT_DATA *b_arry=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*bmpHeight);
	int u,v;
	for(u=0;u<bmpWidth;u++)
	{
		for(v=0;v<bmpHeight;v++)
		{
			b_arry[v].data_R=fp_b_data[v*bmpWidth+u];
			b_arry[v].data_I=0.0;
		}
		FFT_N(b_arry);
		for(v=0;v<bmpHeight;v++)
		{
			fft_b_data[v*bmpWidth+u].data_R=b_arry[v].data_R;
			fft_b_data[v*bmpWidth+u].data_I=b_arry[v].data_I;
		}
	}
	for(v=0;v<bmpHeight;v++)
	{
		FFT_N(fft_b_data+v*bmpWidth);
	}

	//数据变换后存在于R\G\B各自的一维输出数组中
	int r3;
	for(int j=0;j<bmpHeight;j++)
	{
		for(int i=0;i<bmpWidth;i++)
		{
			double temp3=fft_b_data[i+j*bmpWidth].data_R*fft_b_data[i+j*bmpWidth].data_R+fft_b_data[i+j*bmpWidth].data_I*fft_b_data[i+j*bmpWidth].data_I;
			r3=(int)(sqrt(temp3)/100);
			int u1,v1;
			if(i<bmpWidth/2)
				u1=i+bmpWidth/2;
			else
				u1=i-bmpWidth/2;
			if(j<bmpHeight/2)
				v1=j+bmpHeight/2;
			else
				v1=j-bmpHeight/2;
			if(r3>255)
				r3=255;
			if(r3<0)
				r3=0;
			fp_b_data[u1+v1*bmpWidth]=r3;
		}
	}
	//pthread_exit(3);
	pthread_exit(NULL);
}

int main(int argc,char *argv[])
{
	cpu_nub=sysconf(_SC_NPROCESSORS_CONF);	//获得cpu的总数
	printf("总共有%d个CPU\n",cpu_nub);

	FILE *fp_bmp;
	char file_name_input[20];
	double start_time,end_time;
	start_time=clock();

	printf("请输入待处理图片的名字\n");
	gets(file_name_input);				//使用gets()有警告
	fp_bmp=fopen(file_name_input,"rb");
	if(fp_bmp==NULL)
	{
		printf("打开文件失败\n");
		exit(-1);
	}
	unsigned short bfType=0;					//用于判断文件的格式，是否为bmp图像
	fseek(fp_bmp,0L,SEEK_SET);					//相对于开始位置偏移0个字符
	fread(&bfType,sizeof(char),2,fp_bmp);
	if(BM!=bfType)
	{
		printf("打开的不是bmp图像\n");
		exit(-1);
	}

	fseek(fp_bmp,28L,SEEK_SET);
	fread(&bmpBitCount,sizeof(char),2,fp_bmp);		//图像的位深度
	fseek(fp_bmp,10L,SEEK_SET);
	fread(&OffSet,sizeof(char),4,fp_bmp);			//相对图像数据域的偏移长度
	fseek(fp_bmp,18L,SEEK_SET);
	fread(&bmpWidth,sizeof(char),4,fp_bmp);			//图像宽度
	fseek(fp_bmp,22L,SEEK_SET);
	fread(&bmpHeight,sizeof(char),4,fp_bmp);		//图像高度
	printf("位深度=%d\n头部偏移长度=%d\n宽度=%ld\n高度=%ld\n",bmpBitCount,OffSet,bmpWidth,bmpHeight);

	init_temp_array(temp_bmp_arry);					//全局变量－一维数组初始化

	array_line_data=bmpHeight*bmpWidth;				//图像长×宽的值
	Sample_Num=bmpHeight;							//变换点的个数
	pixel_length=line_data*bmpHeight;				//图像包含像素点的总的个数

	/*图像的读入输出*/
	FILE *fp_out=fopen("cpu_fft.bmp","wb+");		//图像的读入输出
	if(fp_out==NULL)
	{
		printf("输出图像打开失败\n");
		exit(-1);
	}
	unsigned char  *fp_data=(unsigned char *)malloc(sizeof(char)*pixel_length);
	unsigned char  *fp_r_data=(unsigned char *)malloc(sizeof(char)*array_line_data);
	unsigned char  *fp_g_data=(unsigned char *)malloc(sizeof(char)*array_line_data);
	unsigned char  *fp_b_data=(unsigned char *)malloc(sizeof(char)*array_line_data);

	fseek(fp_bmp,0L,SEEK_SET);
	fread(fp_data,sizeof(char),OffSet,fp_bmp);
	fwrite(fp_data,sizeof(char),OffSet,fp_out);
	memset(fp_data,0,strlen(fp_data));
	fread(fp_data,sizeof(char),pixel_length,fp_bmp);

	int row,col;
	fseek(fp_bmp,OffSet,SEEK_SET);
	for(row=0;row<bmpHeight;row++)
	{
		fread(temp_bmp_arry,sizeof(char),line_data,fp_bmp);
		for(col=0;col<bmpWidth;col++)
		{
			fp_r_data[row*bmpWidth+col]=temp_bmp_arry[col*3+2];
			fp_g_data[row*bmpWidth+col]=temp_bmp_arry[col*3+1];
			fp_b_data[row*bmpWidth+col]=temp_bmp_arry[col*3+0];
		}
		memset(temp_bmp_arry,0,line_data);
	}

	//首先，把所读的像素数据由实数转化为复数形式，实部置为像素值，虚部置为0

	my_init();			//初始化,第56行
	InitForFFT();

	struct mypara r_para;
	struct mypara g_para;
	struct mypara b_para;

	pthread_t thread1;
	pthread_t thread2;
	pthread_t thread3;

	r_para.thread_id=0;
	r_para.data=fp_r_data;
	g_para.thread_id=1;
	g_para.data=fp_g_data;
	b_para.thread_id=2;
	b_para.data=fp_b_data;

	if(pthread_create(&thread1,NULL,thread_fun_r,&(r_para))!=0)
	{
		perror("pthread_create:");
		exit(-1);
	}

	if(pthread_create(&thread2,NULL,thread_fun_g,&(g_para))!=0)
	{
		perror("pthread_create:");
		exit(-1);
	}

	if(pthread_create(&thread3,NULL,thread_fun_b,&(b_para))!=0)
	{
		perror("pthread_create:");
		exit(-1);
	}

	//线程退出，捕获退出值
	/*
	   int ret1,ret2,ret3;
	   ret1=pthread_join(thread1,(void **)&ret1);
	   ret2=pthread_join(thread2,(void **)&ret2);
	   ret3=pthread_join(thread3,(void **)&ret3);
	   printf("线程一返回值:%d\n",ret1);
	   printf("线程二返回值:%d\n",ret2);
	   printf("线程三返回值:%d\n",ret3);
	   */
	pthread_join(thread1,NULL);
	pthread_join(thread2,NULL);
	pthread_join(thread3,NULL);

	//RGB各自转化后，组成输出数据
	for(row=0;row<bmpHeight;row++)
	{
		memset(temp_bmp_arry,0,line_data);
		fread(temp_bmp_arry,sizeof(char),line_data,fp_bmp);
		for(col=0;col<bmpWidth;col++)
		{
			temp_bmp_arry[col*3+2]=fp_r_data[row*bmpWidth+col];
			temp_bmp_arry[col*3+1]=fp_g_data[row*bmpWidth+col];
			temp_bmp_arry[col*3+0]=fp_b_data[row*bmpWidth+col];
		}
		fwrite(temp_bmp_arry,sizeof(char),line_data,fp_out);
	}

	my_free();
	free(fp_data);
	free(fp_r_data);
	free(fp_g_data);
	free(fp_b_data);

	end_time=clock();
	printf("用时:%lf秒\n",(end_time-start_time)/CLOCKS_PER_SEC);
	return 0;
}
