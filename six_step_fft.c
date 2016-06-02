#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

#define BM 0x4D42					//位图bmp文件的ASCII码
#define PI 3.1415926535898
#define Unit16 unsigned int 

unsigned int OffSet=0;				//位图文件头、信息头的偏移量
long bmpWidth=0;					//图像宽度
long bmpHeight=0;					//图像高度
unsigned int bmpBitCount=0;			//图像的位深度
unsigned char  *temp_bmp_arry;		//存放像素点的一维数组

long int line_data;					//每行应读入的数据的长度
long int array_line_data;			//图像长*宽的长度
long int pixel_length;				//图像所包含的总的像素点的数量

long int Sample_Num_N1;				//转化成矩阵行的长度
long int Sample_Num_N2;				//转化为矩阵列的长度 N1=N2=sqrt(array_line_data)

unsigned char  *fp_data;
unsigned char  *fp_r_data;
unsigned char  *fp_g_data;
unsigned char  *fp_b_data;

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

struct FFT_DATA *fft_r_data;//存放每次变换的数据
struct FFT_DATA *fft_g_data;
struct FFT_DATA *fft_b_data;
struct WN *Wn;
struct FFT_DATA *fft_r_lenth_array;//存放处理结束后的数据
struct FFT_DATA *fft_g_lenth_array;
struct FFT_DATA *fft_b_lenth_array;

struct FFT_DATA Wn_fft;			  //旋转因子结构体——全局变量

struct FFT_DATA **a0_fft_data,**a1_fft_data,**a2_fft_data,**a3_fft_data,**a4_fft_data,
				**a5_fft_data,**c_fft_data;
struct FFT_DATA **g_a0_fft_data,**g_a1_fft_data,**g_a2_fft_data,**g_a3_fft_data,**g_a4_fft_data,
				**g_a5_fft_data,**g_c_fft_data;
struct FFT_DATA **b_a0_fft_data,**b_a1_fft_data,**b_a2_fft_data,**b_a3_fft_data,**b_a4_fft_data,
				**b_a5_fft_data,**b_c_fft_data;

void my_init()
{
	fft_r_data=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N2);//用于存放待变换的数据
	fft_g_data=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N2);
	fft_b_data=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N2);
	Wn=(struct WN *)malloc(sizeof(struct WN)*Sample_Num_N2);
	fft_r_lenth_array=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*array_line_data);
	fft_g_lenth_array=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*array_line_data);
	fft_b_lenth_array=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*array_line_data);

}

void my_free()
{
	free(fft_r_data);
	free(fft_g_data);
	free(fft_b_data);
	free(Wn);
	free(fft_r_lenth_array);
	free(fft_g_lenth_array);
	free(fft_b_lenth_array);
}

void dimen_fft_init()//复数矩阵的初始化
{
	int i;
	a0_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	a1_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	a2_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	a3_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	a4_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	a5_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	c_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	g_a0_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	g_a1_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	g_a2_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	g_a3_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	g_a4_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	g_a5_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	g_c_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	b_a0_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	b_a1_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	b_a2_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	b_a3_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	b_a4_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	b_a5_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	b_c_fft_data=(struct FFT_DATA **)malloc(sizeof(struct FFT_DATA *)*Sample_Num_N2);
	for(i=0;i<Sample_Num_N2;i++)
	{
		a0_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		a1_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		a2_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		a3_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		a4_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		a5_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		c_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		g_a0_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		g_a1_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		g_a2_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		g_a3_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		g_a4_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		g_a5_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		g_c_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		b_a0_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		b_a1_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		b_a2_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		b_a3_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		b_a4_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		b_a5_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
		b_c_fft_data[i]=(struct FFT_DATA *)malloc(sizeof(struct FFT_DATA)*Sample_Num_N1);
	}
}

void dimen_fft_free()//复数矩阵的释放
{
	int i;
	for(i=0;i<Sample_Num_N2;i++)
	{
		free(a0_fft_data[i]);
		free(a1_fft_data[i]);
		free(a2_fft_data[i]);
		free(a3_fft_data[i]);
		free(a4_fft_data[i]);
		free(a5_fft_data[i]);
		free(c_fft_data[i]);
		free(g_a0_fft_data[i]);
		free(g_a1_fft_data[i]);
		free(g_a2_fft_data[i]);
		free(g_a3_fft_data[i]);
		free(g_a4_fft_data[i]);
		free(g_a5_fft_data[i]);
		free(g_c_fft_data[i]);
		free(b_a0_fft_data[i]);
		free(b_a1_fft_data[i]);
		free(b_a2_fft_data[i]);
		free(b_a3_fft_data[i]);
		free(b_a4_fft_data[i]);
		free(b_a5_fft_data[i]);
		free(b_c_fft_data[i]);
	}
	free(a0_fft_data);
	free(a1_fft_data);
	free(a2_fft_data);
	free(a3_fft_data);
	free(a4_fft_data);
	free(a5_fft_data);
	free(c_fft_data);
	free(g_a0_fft_data);
	free(g_a1_fft_data);
	free(g_a2_fft_data);
	free(g_a3_fft_data);
	free(g_a4_fft_data);
	free(g_a5_fft_data);
	free(g_c_fft_data);
	free(b_a0_fft_data);
	free(b_a1_fft_data);
	free(b_a2_fft_data);
	free(b_a3_fft_data);
	free(b_a4_fft_data);
	free(b_a5_fft_data);
	free(b_c_fft_data);
}

void array_to_complex_array(unsigned char *fp_rgb_data,struct FFT_DATA **a0)//把一维实数矩阵形式，转化为二维复数矩阵形式
{
	int i,j,k;
	k=0;
	for(j=0;j<Sample_Num_N1;j++)
	{
		for(i=0;i<Sample_Num_N2;i++)
		{
			a0[i][j].data_R=fp_rgb_data[k];
			a0[i][j].data_I=0.0;
			k=k+1;
		}
	}
}

void zero_complex_array(struct FFT_DATA *arry)		//一维复数数组置为零
{
	int i;
	for(i=0;i<Sample_Num_N2;i++)
	{
		arry[i].data_R=0.0;
		arry[i].data_I=0.0;
	}
}

void  complex_array_transpose(struct FFT_DATA **array1,struct FFT_DATA **array2)//复数矩阵的转置
{
	int i,j;
	for(i=0;i<Sample_Num_N2;i++)
	{
		for(j=0;j<Sample_Num_N1;j++)
		{
			array2[j][i].data_R=array1[i][j].data_R;
			array2[j][i].data_I=array1[i][j].data_I;
		}
	}
}


void InitForFFT(long int N)						
{
	int i;
	for(i=0;i<N;i++)			//(2*PI*n/N)——与公式相比，缺少一个k
	{
		Wn[i].sin=sin(PI*2*i/N);
		Wn[i].cos=cos(PI*2*i/N);
	}
}

void  twist_coefficient(struct FFT_DATA temp_fft,float tp,int t)	//实现复数的连乘
{	
	struct FFT_DATA my_temp_fft;
	for(;t<Sample_Num_N2;t++)
	{
		my_temp_fft.data_R=temp_fft.data_R*cos(tp*t)+temp_fft.data_I*sin(tp*t);
		my_temp_fft.data_I=temp_fft.data_I*cos(tp*t)-temp_fft.data_R*sin(tp*t);
		temp_fft.data_R=my_temp_fft.data_R;
		temp_fft.data_I=my_temp_fft.data_I;
	}
	Wn_fft.data_R=temp_fft.data_R;
	Wn_fft.data_I=temp_fft.data_I;
}

void FFT_N1(struct FFT_DATA *fft_data,long int number,int col,struct FFT_DATA **a2,struct FFT_DATA **a3) //R G B分别进行FFT变换
{
	Unit16 xx,Num;
	Unit16 i,j,k,b,p,L;
	float TR,TI,temp;
	Num=(Unit16)(log(number)/log(2)+0.0001);	//num=log(number)
	//数据的倒序与原位
	for(i=0;i<number;i++)
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

	for(i=0;i<number;i++)
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
			for(k=j;k<number;k=k+2*b)//计算每一级蝶形运算的旋转因子,第k与第k+b距离的两个点进行蝶形											运算
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
	//变换完成后的数据存放在一维数组fft_data中
	for(int t=0;t<number;t++)
	{
		a2[t][col].data_R=fft_data[t].data_R;
		a2[t][col].data_I=fft_data[t].data_I;
	}
	//乘以一个旋转因子
	struct FFT_DATA temp_fft_data;	//辅助复数变量
	for(int k2=0;k2<number;k2++)
	{
		float temp=(2*PI*k2)/array_line_data;
		temp_fft_data.data_R=cos(temp);
		temp_fft_data.data_I=sin(temp);
		
		Wn_fft.data_R=0.0;
		Wn_fft.data_I=0.0;

		twist_coefficient(temp_fft_data,temp,2);//求得旋转因袭,保存在Wn_fft中

		//a2与旋转因子Wn_fft相乘
		a3[k2][col].data_R=a2[k2][col].data_R*Wn_fft.data_R+a2[k2][col].data_I*Wn_fft.data_I;		
		a3[k2][col].data_I=a2[k2][col].data_I*Wn_fft.data_R-a2[k2][col].data_R*Wn_fft.data_I;
	}
}


void FFT_N2(struct FFT_DATA *fft_data,long int number,int col,struct FFT_DATA **a5) //R G B分别进行FFT变换
{
	Unit16 xx,Num;
	Unit16 i,j,k,b,p,L;
	float TR,TI,temp;
	Num=(Unit16)(log(number)/log(2)+0.0001);	//num=log(number)

	//数据的倒序与原位
	for(i=0;i<number;i++)
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

	for(i=0;i<number;i++)
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
			for(k=j;k<number;k=k+2*b)//计算每一级蝶形运算的旋转因子,第k与第k+b距离的两个点进行蝶形											运算
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
	//变换完成后的数据存放在一维数组fft_data中
	for(int t=0;t<number;t++)
	{
		a5[t][col].data_R=fft_data[t].data_R;
		a5[t][col].data_I=fft_data[t].data_I;
	}
	
}


int main(int argc,char *argv[])
{
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
	pixel_length=line_data*bmpHeight;				//图像包含像素点的总的个数

	Sample_Num_N1=(Unit16)(sqrt(array_line_data));	//复数矩阵行的长度
	Sample_Num_N2=(Unit16)(sqrt(array_line_data)); //复数矩阵列的长度

	printf("N1=%ld\tN2=%ld\n",Sample_Num_N1,Sample_Num_N2);
	/*图像的读入输出*/
	FILE *fp_out=fopen("six_lena.bmp","wb+");		//图像的读入输出
	if(fp_out==NULL)
	{
		printf("输出图像打开失败\n");
		exit(-1);
	}

	fp_data=(unsigned char *)malloc(sizeof(char)*pixel_length);
	fp_r_data=(unsigned char *)malloc(sizeof(char)*array_line_data);
	fp_g_data=(unsigned char *)malloc(sizeof(char)*array_line_data);
	fp_b_data=(unsigned char *)malloc(sizeof(char)*array_line_data);

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

	dimen_fft_init();	//复数矩阵的初始化
	my_init();			//复数一维数组的初始化

	//六步FFT变换
	//六步快速傅里叶变换——数据准备阶段，由一维实数数据转变成二维复数矩阵形式
	array_to_complex_array(fp_r_data,a0_fft_data);		//形参R、G、B一维数组形式的数据
	//第一步：矩阵转置
	complex_array_transpose(a0_fft_data,a1_fft_data);//转置后的数据存放在a1_fft_data[][]数组中
	//第二步：第一次傅里叶变换
	InitForFFT(Sample_Num_N2);			//sin与cos的初始化
	for(int j1=0;j1<Sample_Num_N1;j1++)
	{
		for(int j2=0;j2<Sample_Num_N2;j2++)
		{
			fft_r_data[j2].data_R=a1_fft_data[j2][j1].data_R;	//取一列数据
			fft_r_data[j2].data_I=a1_fft_data[j2][j1].data_I;
		}
		FFT_N1(fft_r_data,Sample_Num_N2,j1,a2_fft_data,a3_fft_data);	//一列数据进行变换
		zero_complex_array(fft_r_data);
	}

	printf("红色信息经过傅里叶变换后，数据是:\n");
	for(int i=0;i<5;i++)
	{
		for(int j=0;j<3;j++)
		{
			printf("%f\t%tf\t",a3_fft_data[i][j].data_R,a3_fft_data[i][j].data_I);
		}
		printf("\n");
	}

	//第四步：转置
	complex_array_transpose(a3_fft_data,a4_fft_data);
	//第五步：第二次进行快速傅里叶变换
	zero_complex_array(fft_r_data);
	for(int k2=0;k2<Sample_Num_N2;k2++)
	{
		for(int j1=0;j1<Sample_Num_N1;j1++)
		{
			fft_r_data[j1].data_R=a4_fft_data[j1][k2].data_R;//取一列数据
			fft_r_data[j1].data_I=a4_fft_data[j1][k2].data_I;
		}
		FFT_N2(fft_r_data,Sample_Num_N1,k2,a5_fft_data);
		zero_complex_array(fft_r_data);
	}

	//第六步：数组转置
	complex_array_transpose(a5_fft_data,c_fft_data);//变化完成结果存放在复数矩阵c_fft_data中
	//数据后处理：把处理结束后，复数矩阵中的数据转化成复数一维数组中
	int r_k=0;
	for(int i=0;i<Sample_Num_N2;i++)
	{
		for(int j=0;j<Sample_Num_N1;j++)
		{	
			fft_r_lenth_array[r_k].data_R=c_fft_data[j][i].data_R;
			fft_r_lenth_array[r_k].data_I=c_fft_data[j][i].data_I;
			r_k=r_k+1;
		}
	}		//红色信息数据处理完毕，处理后数据存放在fft_r_lenth_array一维数组中
	printf("红色数据处理完成\n");
	printf("处理完成，红色后10个数据是:\n");
	for(int i=512*512;i>(512*512-10);i--)
	{
		printf("%f\t%f\n",fft_r_lenth_array[i].data_R,fft_r_lenth_array[i].data_I);
	}

	//六步快速傅里叶变换——数据准备阶段，由一维实数数据转变成二维复数矩阵形式
	array_to_complex_array(fp_g_data,g_a0_fft_data);		//形参R、G、B一维数组形式的数据
	//第一步：矩阵转置
	complex_array_transpose(g_a0_fft_data,g_a1_fft_data);//转置后的数据存放在a1_fft_data[][]数组中
	//第二步：第一次傅里叶变换
	InitForFFT(Sample_Num_N2);			//sin与cos的初始化
	for(int j1=0;j1<Sample_Num_N1;j1++)
	{
		for(int j2=0;j2<Sample_Num_N2;j2++)
		{
			fft_g_data[j2].data_R=g_a1_fft_data[j2][j1].data_R;	//取一列数据
			fft_g_data[j2].data_I=g_a1_fft_data[j2][j1].data_I;
		}
		FFT_N1(fft_g_data,Sample_Num_N2,j1,g_a2_fft_data,g_a3_fft_data);	//一列数据进行变换
		zero_complex_array(fft_g_data);
	}

	printf("绿色信息经过傅里叶变换后，数据是:\n");
	for(int i=0;i<5;i++)
	{
		for(int j=0;j<3;j++)
		{
			printf("%f\t%tf\t",g_a3_fft_data[i][j].data_R,g_a3_fft_data[i][j].data_I);
		}
		printf("\n");
	}

	//第四步：转置
	complex_array_transpose(g_a3_fft_data,g_a4_fft_data);
	//第五步：第二次进行快速傅里叶变换
	zero_complex_array(fft_g_data);
	for(int k2=0;k2<Sample_Num_N2;k2++)
	{
		for(int j1=0;j1<Sample_Num_N1;j1++)
		{
			fft_g_data[j1].data_R=g_a4_fft_data[j1][k2].data_R;//取一列数据
			fft_g_data[j1].data_I=g_a4_fft_data[j1][k2].data_I;
		}
		FFT_N2(fft_g_data,Sample_Num_N1,k2,g_a5_fft_data);
		zero_complex_array(fft_g_data);
	}

	//第六步：数组转置
	complex_array_transpose(g_a5_fft_data,g_c_fft_data);//变化完成结果存放在复数矩阵c_fft_data中
	//数据后处理：把处理结束后，复数矩阵中的数据转化成复数一维数组中
	int g_k=0;
	for(int i=0;i<Sample_Num_N2;i++)
	{
		for(int j=0;j<Sample_Num_N1;j++)
		{	
			fft_g_lenth_array[g_k].data_R=g_c_fft_data[j][i].data_R;
			fft_g_lenth_array[g_k].data_I=g_c_fft_data[j][i].data_I;
			g_k=g_k+1;
		}
	}		//绿色信息数据处理完毕，处理后数据存放在fft_g_lenth_array一维数组中
	printf("绿色数据处理完成\n");
	printf("处理完成，绿色后10个数据是:\n");
	for(int i=512*512;i>(512*512-10);i--)
	{
		printf("%f\t%f\n",fft_g_lenth_array[i].data_R,fft_g_lenth_array[i].data_I);
	}

	//六步快速傅里叶变换——数据准备阶段，由一维实数数据转变成二维复数矩阵形式
	array_to_complex_array(fp_b_data,b_a0_fft_data);		//形参R、G、B一维数组形式的数据
	//第一步：矩阵转置
	complex_array_transpose(b_a0_fft_data,b_a1_fft_data);//转置后的数据存放在a1_fft_data[][]数组中
	//第二步：第一次傅里叶变换
	InitForFFT(Sample_Num_N2);			//sin与cos的初始化
	for(int j1=0;j1<Sample_Num_N1;j1++)
	{
		for(int j2=0;j2<Sample_Num_N2;j2++)
		{
			fft_b_data[j2].data_R=b_a1_fft_data[j2][j1].data_R;	//取一列数据
			fft_b_data[j2].data_I=b_a1_fft_data[j2][j1].data_I;
		}
		FFT_N1(fft_b_data,Sample_Num_N2,j1,b_a2_fft_data,b_a3_fft_data);	//一列数据进行变换
		zero_complex_array(fft_b_data);
	}

	printf("蓝色信息经过傅里叶变换后，数据是:\n");
	for(int i=0;i<5;i++)
	{
		for(int j=0;j<3;j++)
		{
			printf("%f\t%tf\t",b_a3_fft_data[i][j].data_R,b_a3_fft_data[i][j].data_I);
		}
		printf("\n");
	}

	//第四步：转置
	complex_array_transpose(b_a3_fft_data,b_a4_fft_data);
	//第五步：第二次进行快速傅里叶变换
	zero_complex_array(fft_b_data);
	for(int k2=0;k2<Sample_Num_N2;k2++)
	{
		for(int j1=0;j1<Sample_Num_N1;j1++)
		{
			fft_b_data[j1].data_R=b_a4_fft_data[j1][k2].data_R;//取一列数据
			fft_b_data[j1].data_I=b_a4_fft_data[j1][k2].data_I;
		}
		FFT_N2(fft_b_data,Sample_Num_N1,k2,b_a5_fft_data);
		zero_complex_array(fft_b_data);
	}

	//第六步：数组转置
	complex_array_transpose(b_a5_fft_data,b_c_fft_data);//变化完成结果存放在复数矩阵c_fft_data中
	//数据后处理：把处理结束后，复数矩阵中的数据转化成复数一维数组中
	int b_k=0;
	for(int i=0;i<Sample_Num_N2;i++)
	{
		for(int j=0;j<Sample_Num_N1;j++)
		{	
			fft_b_lenth_array[b_k].data_R=b_c_fft_data[j][i].data_R;
			fft_b_lenth_array[b_k].data_I=b_c_fft_data[j][i].data_I;
			b_k=b_k+1;
		}
	}		//红色信息数据处理完毕，处理后数据存放在fft_r_lenth_array一维数组中
	printf("蓝色数据处理完成\n");
	printf("处理完成，蓝色后10个数据是:\n");
	for(int i=512*512;i>(512*512-10);i--)
	{
		printf("%f\t%f\n",fft_b_lenth_array[i].data_R,fft_b_lenth_array[i].data_I);
	}

	//数据在傅里叶变换后存在于R、G、B各自的一维输出数组中
	int r1,r2,r3;
	for(int j=0;j<bmpHeight;j++)
	{
		for(int i=0;i<bmpWidth;i++)
		{
			double temp1=fft_r_lenth_array[i+j*bmpWidth].data_R*fft_r_lenth_array[i+j*bmpWidth].data_R+fft_r_lenth_array[i+j*bmpWidth].data_I*fft_r_lenth_array[i+j*bmpWidth].data_I;

			double temp2=fft_g_lenth_array[i+j*bmpWidth].data_R*fft_g_lenth_array[i+j*bmpWidth].data_R+fft_g_lenth_array[i+j*bmpWidth].data_I*fft_g_lenth_array[i+j*bmpWidth].data_I;

			double temp3=fft_b_lenth_array[i+j*bmpWidth].data_R*fft_b_lenth_array[i+j*bmpWidth].data_R+fft_b_lenth_array[i+j*bmpWidth].data_I*fft_b_lenth_array[i+j*bmpWidth].data_I;

			r1=(int)(sqrt(temp1)/100);
			r2=(int)(sqrt(temp2)/100);
			r3=(int)(sqrt(temp3)/100);
			int u,v;
			if(i<bmpWidth/2)
				u=i+bmpWidth/2;
			else
				u=i-bmpWidth/2;
			if(j<bmpHeight/2)
				v=j+bmpHeight/2;
			else
				v=j-bmpHeight/2;

			if(r1>255)
				r1=255;
			if(r1<0)
				r1=0;
			fp_r_data[u+v*bmpWidth]=r1;

			if(r2>255)
				r2=255;
			if(r2<0)			
				r2=0;
			fp_g_data[u+v*bmpWidth]=r2;

			if(r3>255)
				r3=255;
			if(r3<0)
				r3=0;
			fp_b_data[u+v*bmpWidth]=r3;
		}
	}

	printf("RGB数据变换完成，前10个数据是\n");
	for(int i=0;i<10;i++)
	{
		printf("%d\t%d\t%d\n",fp_r_data[i],fp_g_data[i],fp_b_data[i]);
	}
	printf("最后10个数据是\n");
	for(int i=512*512;i>(512*512-10);i--)
	{
		printf("%d\t%d\t%d\n",fp_r_data[i],fp_g_data[i],fp_b_data[i]);
	}

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

	//fwrite(fp_data,sizeof(char),pixel_length,fp_out);//输出为原图——用于测试
	my_free();
	dimen_fft_free();
	free(fp_data);
	free(fp_r_data);
	free(fp_g_data);
	free(fp_b_data);
	
	end_time=clock();
	printf("用时:%lf秒\n",(end_time-start_time)/CLOCKS_PER_SEC);
	return 0;
}
