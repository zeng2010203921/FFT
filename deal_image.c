#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

#define bsize   8   //DCT变换块的大小　8*8
#define BM 0x4D42 //位图bmp文件的ASCII码
#define PI 3.1415926535898

unsigned int OffSet=0;//位图文件头、信息头的偏移量
long bmpWidth=0;//图像宽度
long bmpHeight=0;//图像高度
unsigned int bmpBitCount=0;//图像的位深度

long int line_data;	//实际每行数据的长度

double  **r_array;  //存放rgb数组
double  **g_array;
double  **b_array;
unsigned char  *temp_bmp_arry;

unsigned char  **laplase_r;
unsigned char  **laplase_g;
unsigned char  **laplase_b;

double **dct_r_array;//存放dct变换的数组
double **dct_g_array;
double **dct_b_array;

double **idct_r_array;//存放idct变换数组
double **idct_g_array;
double **idct_b_array;

void my_malloc()//自定义动态分配二位数组
{
	int i;
	r_array=(double **)malloc(sizeof(double *)*bmpHeight);
	g_array=(double **)malloc(sizeof(double *)*bmpHeight);
	b_array=(double **)malloc(sizeof(double *)*bmpHeight);

	laplase_r=(unsigned char **)malloc(sizeof(unsigned char*)*bmpHeight);
	laplase_g=(unsigned char **)malloc(sizeof(unsigned char*)*bmpHeight);
	laplase_b=(unsigned char **)malloc(sizeof(unsigned char*)*bmpHeight);
	
	dct_r_array=(double **)malloc(sizeof(double *)*bmpHeight);
	dct_g_array=(double **)malloc(sizeof(double *)*bmpHeight);
	dct_b_array=(double **)malloc(sizeof(double *)*bmpHeight);

	idct_r_array=(double **)malloc(sizeof(double *)*bmpHeight);
	idct_g_array=(double **)malloc(sizeof(double *)*bmpHeight);
	idct_b_array=(double **)malloc(sizeof(double *)*bmpHeight);

	for(i=0;i<bmpHeight;i++)
	{
		r_array[i]=(double *)malloc(sizeof(double)*bmpWidth);
		g_array[i]=(double *)malloc(sizeof(double)*bmpWidth);
		b_array[i]=(double *)malloc(sizeof(double)*bmpWidth);
	
		laplase_r[i]=(unsigned char *)malloc(sizeof(unsigned char)*bmpWidth);
		laplase_g[i]=(unsigned char *)malloc(sizeof(unsigned char)*bmpWidth);
		laplase_b[i]=(unsigned char *)malloc(sizeof(unsigned char)*bmpWidth);

		dct_r_array[i]=(double *)malloc(sizeof(double)*bmpWidth);
		dct_g_array[i]=(double *)malloc(sizeof(double)*bmpWidth);
		dct_b_array[i]=(double *)malloc(sizeof(double)*bmpWidth);
		
		idct_r_array[i]=(double *)malloc(sizeof(double)*bmpWidth);
		idct_g_array[i]=(double *)malloc(sizeof(double)*bmpWidth);
		idct_b_array[i]=(double *)malloc(sizeof(double)*bmpWidth);
	}
}

void my_free()
{
	int i;
	for(i=0;i<bmpHeight;i++)
	{
		free(r_array[i]);
		free(g_array[i]);
		free(b_array[i]);
		free(laplase_r[i]);
		free(laplase_g[i]);
		free(laplase_b[i]);
		free(dct_r_array[i]);
		free(dct_g_array[i]);
		free(dct_b_array[i]);
		free(idct_r_array[i]);
		free(idct_g_array[i]);
		free(idct_b_array[i]);
	}
	free(r_array);
	free(g_array);
	free(b_array);
	free(laplase_r);
	free(laplase_g);
	free(laplase_b);
	free(dct_r_array);
	free(dct_g_array);
	free(dct_b_array);
	free(idct_r_array);
	free(idct_g_array);
	free(idct_b_array);
}

void init_temp_array(double *p)//对中间数组进行初始化
{
	line_data=(bmpWidth*bmpBitCount+31)/8;
	line_data=line_data/4*4;
	temp_bmp_arry=(unsigned char *)malloc(sizeof(unsigned char)*line_data);
}

void coeff(double **dct_coef,int n)   //自定义函数，初始化DCT系数
{
	double sqrt_1=1.0/sqrt(2.0);
	int i,j;
	for(i=0;i<n;i++)
	{
		dct_coef[0][i]=sqrt_1;
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			dct_coef[i][j]=cos(i*PI*(j+0.5)/(double)n);
		}
	}
}

void dct(double **a,double **b,double **c,int n)
{
	double x;
	double af[n][n];
	int i,j,k;
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			x=0.0;
			for(k=0;k<n;k++)
			{
				x+=a[i][k]*b[k][j];
			}
			af[i][j]=x;
		}
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			x=0.0;
			for(k=0;k<n;k++)
			{
				x+=c[i][k]*af[k][j];
			}
			a[i][j]=2.0*x/((double)n);
		}
	}
}

/*图像进行DCT变换*/

void dctTrans(double **img,double **dct_image,long int iw,long int ih,int type)
{
	int i,j;
	int k,l;
	int iter_num=bmpHeight/bsize;
	
	double **dct_coef=(double **)malloc(sizeof(double *)*bsize);
	double **dct_coeft=(double **)malloc(sizeof(double *)*bsize);
	double **image=(double **)malloc(sizeof(double *)*bsize);
	for(i=0;i<bsize;i++)
	{
		dct_coef[i]=(double *)malloc(sizeof(double)*bsize);
		dct_coeft[i]=(double *)malloc(sizeof(double)*bsize);
		image[i]=(double *)malloc(sizeof(double)*bsize);
	}

	coeff(dct_coef,bsize);//自定义函数，初始化dct系数
	//定义转置矩阵系数——逆DCT变换系数
	
	for(i=0;i<bsize;i++)
	{
		for(j=0;j<bsize;j++)
		{
			dct_coeft[i][j]=dct_coef[j][i];
		}
	}
	
	if(type==1)//type==1进行DCT变换
	{
		for(i=0;i<iter_num;i++)
		{
			for(j=0;j<iter_num;j++)
			{	
       	    			//取bsize*bsize的图像块，存放在image[][]
				for(k=0;k<bsize;k++)
				{
   	         				for(l=0;l<bsize;l++)
					{
						image[k][l]=img[i*bsize+k][j*bsize+l];//把整副图像的某一块提取出来
					}
				}
				
				//对bsize*bsize块进行DCT变换
				dct(image,dct_coeft,dct_coef,bsize);
				
				//输出DCT变换后的图像
				for(k=0;k<bsize;k++)
				{
 	          				for(l=0;l<bsize;l++)
					{
						dct_image[i*bsize+k][j*bsize+l]=image[k][l];
					}
				}
			}
		}
	}
	else	  //进行逆DCT变换
	{
		
		for(i=0;i<iter_num;i++)
		{
			for(j=0;j<iter_num;j++)
			{	
   	        			//取bsize*bsize的图像块，存放在image[][]
				for(k=0;k<bsize;k++)
				{
 	            	for(l=0;l<bsize;l++)
					{
						image[k][l]=img[i*bsize+k][j*bsize+l];
					}
				}
				
				//对bsize*bsize块进行逆DCT变换
				dct(image,dct_coef,dct_coeft,bsize);
				
				//输出DCT变换后的图像
				for(k=0;k<bsize;k++)
				{
 	           		for(l=0;l<bsize;l++)
					{
						dct_image[i*bsize+k][j*bsize+l]=image[k][l];
					}
				}
			}
		}
	}
	for(i=0;i<bsize;i++)
	{
		free(dct_coef[i]);
		free(dct_coeft[i]);
		free(image[i]);
	}
	free(dct_coef);
	free(dct_coeft);
	free(image);

	//return dct_image;

}


int main(int argc,char *argv[])
{
	FILE *fp_bmp;
	char file_name_input[20];
	double start_time,end_time;
	start_time=clock();

	printf("请输入待处理图片的名字\n");
	gets(file_name_input);//使用gets()有警告
	fp_bmp=fopen(file_name_input,"rb");
	if(fp_bmp==NULL)
	{
		printf("打开文件失败\n");
		exit(-1);
	}

	unsigned short bfType=0;//用于判断文件的格式，是否为bmp图像
	fseek(fp_bmp,0L,SEEK_SET);//相对于开始位置偏移0个字符
	fread(&bfType,sizeof(char),2,fp_bmp);
	if(BM!=bfType)
	{
		printf("打开的不是bmp图像\n");
		exit(-1);
	}
	
	fseek(fp_bmp,28L,SEEK_SET);
	fread(&bmpBitCount,sizeof(char),2,fp_bmp);//图像的位深度
	fseek(fp_bmp,10L,SEEK_SET);
	fread(&OffSet,sizeof(char),4,fp_bmp);//相对图像数据域的偏移长度
	fseek(fp_bmp,18L,SEEK_SET);
	fread(&bmpWidth,sizeof(char),4,fp_bmp);//图像宽度
	fseek(fp_bmp,22L,SEEK_SET);
	fread(&bmpHeight,sizeof(char),4,fp_bmp);//图像高度
	printf("位深度=%d\n头部偏移长度=%d\n宽度=%ld\n高度=%ld\n",bmpBitCount,OffSet,bmpWidth,bmpHeight);
	
	init_temp_array(temp_bmp_arry);//全局变量－一维数组初始化
	/*图像的读入输出*/
	FILE *fp_out=fopen("lena_pixel.bmp","wb+");//图像的读入输出
	if(fp_out==NULL)
	{
		printf("输出图像打开失败\n");
		exit(-1);
	}
	long int pixel_length=line_data*bmpHeight;
	unsigned char  *fp_data=(unsigned char *)malloc(sizeof(char)*pixel_length);
	printf("数组长度:%d\n",pixel_length);
	fseek(fp_bmp,0L,SEEK_SET);
	fread(fp_data,sizeof(char),OffSet,fp_bmp);
	fwrite(fp_data,sizeof(char),OffSet,fp_out);
	memset(fp_data,0,strlen(fp_data));

	fseek(fp_bmp,OffSet,SEEK_SET);
	fread(fp_data,sizeof(char),pixel_length,fp_bmp);
	srand(time(NULL));
	for(int i=0;i<pixel_length;i++)
	{
		fp_data[i]=fp_data[i]-rand()%128;
		if(fp_data[i]<0)
			fp_data[i]=0;
	}
	fwrite(fp_data,sizeof(char),pixel_length,fp_out);

	/*灰度图像的输出*/
	FILE *fp_gray_out=fopen("lena_gray.bmp","wb+");//灰度图像的输出
	if(fp_gray_out==NULL)
	{
		printf("灰度图像输出文件打开失败\n");
		exit(-1);
	}
	long int gray_lenth=bmpWidth*line_data;
	unsigned char  *fp_gray_data=(unsigned char *)malloc(sizeof(char)*gray_lenth);//灰度图像处理
	fseek(fp_bmp,0L,SEEK_SET);
	fread(fp_gray_data,sizeof(char),OffSet,fp_bmp);
	fwrite(fp_gray_data,sizeof(char),OffSet,fp_gray_out);
	memset(fp_gray_data,0,strlen(fp_gray_data));

	int row,col;
	for(row=0;row<bmpHeight*3;row++)//进行灰度变换
	{
		for(col=0;col<bmpWidth;col++)
		{
			fp_gray_data[row*bmpWidth+col+2]=(unsigned char)(0.3*fp_data[row*bmpWidth+col+2]+0.59*fp_data[row*bmpWidth+col+1]+0.11*fp_data[row*bmpWidth+col]);//B

			fp_gray_data[row*bmpWidth+col+1]=(unsigned char)(0.3*fp_data[row*bmpWidth+col+2]+0.59*fp_data[row*bmpWidth+col+1]+0.11*fp_data[row*bmpWidth+col]);//G

			fp_gray_data[row*bmpWidth+col+0]=(unsigned char)(0.3*fp_data[row*bmpWidth+col+2]+0.59*fp_data[row*bmpWidth+col+1]+0.11*fp_data[row*bmpWidth+col]);//R
		}
	}
	fseek(fp_gray_out,OffSet,SEEK_SET);
	fwrite(fp_gray_data,sizeof(char),gray_lenth,fp_gray_out);

	my_malloc();//全局变量－二维数组初始化

	fseek(fp_bmp,OffSet,SEEK_SET);
	//读数据到r_array[][],g_array[][],b_array[][]中
	for(row=0;row<bmpHeight;row++)
	{
		fread(temp_bmp_arry,sizeof(char),line_data,fp_bmp);
		for(col=0;col<bmpWidth;col++)
		{
			r_array[bmpHeight-1-row][col]=temp_bmp_arry[col*3+2];
			g_array[bmpHeight-1-row][col]=temp_bmp_arry[col*3+1];
			b_array[bmpHeight-1-row][col]=temp_bmp_arry[col*3+0];
		}
		memset(temp_bmp_arry,0,line_data);
	}

	/*拉普拉斯图像输出*/
	FILE *fp_laplase_out=fopen("lena_laplase.bmp","wb+");//拉普拉斯图像的输出
	if(fp_laplase_out==NULL)
	{
		printf("拉普拉斯图像输出文件打开失败\n");
		exit(-1);
	}
	unsigned char *laplase_data=(unsigned char *)malloc(sizeof(char)*OffSet);
	fseek(fp_bmp,0L,SEEK_SET);
	fread(laplase_data,sizeof(char),OffSet,fp_bmp);
	fwrite(laplase_data,sizeof(char),OffSet,fp_laplase_out);
	free(laplase_data);
	
	//第一行、最后一行不做处理
	for(col=0;col<bmpWidth;col++)
	{
		laplase_r[0][col]=(unsigned char)r_array[0][col];
		laplase_g[0][col]=(unsigned char)g_array[0][col];
		laplase_b[0][col]=(unsigned char)b_array[0][col];
		laplase_r[bmpHeight-1][col]=(unsigned char)r_array[bmpHeight-1][col];
		laplase_g[bmpHeight-1][col]=(unsigned char)g_array[bmpHeight-1][col];
		laplase_b[bmpHeight-1][col]=(unsigned char)b_array[bmpHeight-1][col];
	}
	//第一列、最后一列不做处理
	for(row=1;row<bmpHeight-1;row++)
	{
		laplase_r[row][0]=(unsigned char)r_array[row][0];
		laplase_g[row][0]=(unsigned char)g_array[row][0];
		laplase_b[row][0]=(unsigned char)b_array[row][0];
		laplase_r[row][bmpWidth-1]=(unsigned char)r_array[row][bmpWidth-1];
		laplase_g[row][bmpWidth-1]=(unsigned char)g_array[row][bmpWidth-1];
		laplase_b[row][bmpWidth-1]=(unsigned char)b_array[row][bmpWidth-1];
	}
	//其他数据块
	for(row=1;row<bmpHeight-1;row++)
	{
		for(col=1;col<bmpWidth-1;col++)
		{
			laplase_r[row][col]=(unsigned char)(r_array[row+1][col]+r_array[row-1][col]+r_array[row][col+1]+r_array[row][col-1]-4*r_array[row][col]);
			laplase_g[row][col]=(unsigned char)(g_array[row+1][col]+g_array[row-1][col]+g_array[row][col+1]+g_array[row][col-1]-4*g_array[row][col]);
			laplase_b[row][col]=(unsigned char)(b_array[row+1][col]+b_array[row-1][col]+b_array[row][col+1]+b_array[row][col-1]-4*b_array[row][col]);
		}
	}
	unsigned char  *laplase_pixel_out=(unsigned char *)malloc(sizeof(char)*line_data);
	for(row=0;row<bmpHeight;row++)
	{
		for(col=0;col<bmpWidth;col++)
		{
			laplase_pixel_out[col*3+2]=(unsigned char)laplase_r[bmpHeight-1-row][col];
			laplase_pixel_out[col*3+1]=(unsigned char)laplase_g[bmpHeight-1-row][col];
			laplase_pixel_out[col*3+0]=(unsigned char)laplase_b[bmpHeight-1-row][col];
		}
		fwrite(laplase_pixel_out,sizeof(char),line_data,fp_laplase_out);
	}


	/*图像进行DCT－离散余弦变换*/
	
	dctTrans(r_array,dct_r_array,bmpWidth,bmpHeight,1);
	dctTrans(g_array,dct_g_array,bmpWidth,bmpHeight,1);
	dctTrans(b_array,dct_b_array,bmpWidth,bmpHeight,1);
	
	FILE *fp_dct_out=fopen("lena_dct.bmp","wb+");//拉普拉斯图像的输出
	if(fp_dct_out==NULL)
	{
		printf("离散余弦变换输出文件打开失败\n");
		exit(-1);
	}
	unsigned char *dct_data=(unsigned char *)malloc(sizeof(char)*OffSet);
	fseek(fp_bmp,0L,SEEK_SET);
	fread(dct_data,sizeof(char),OffSet,fp_bmp);
	fwrite(dct_data,sizeof(char),OffSet,fp_dct_out);
	free(dct_data);
	unsigned char  *dct_pixel_out=(unsigned char *)malloc(sizeof(char)*line_data);
	for(row=0;row<bmpHeight;row++)
	{
		for(col=0;col<bmpWidth;col++)
		{
			dct_pixel_out[col*3+2]=(unsigned char)dct_r_array[bmpHeight-1-row][col];
			dct_pixel_out[col*3+1]=(unsigned char)dct_g_array[bmpHeight-1-row][col];
			dct_pixel_out[col*3+0]=(unsigned char)dct_b_array[bmpHeight-1-row][col];
		}
		fwrite(dct_pixel_out,sizeof(char),line_data,fp_dct_out);
	}	
	
	
	//逆DCT变换,有一个量化步数需要确定
	dctTrans(dct_r_array,idct_r_array,bmpWidth,bmpHeight,2);
	dctTrans(dct_g_array,idct_g_array,bmpWidth,bmpHeight,2);
	dctTrans(dct_b_array,idct_b_array,bmpWidth,bmpHeight,2);
	
	FILE *fp_idct_out=fopen("lena_idct.bmp","wb+");//拉普拉斯图像的输出
	if(fp_idct_out==NULL)
	{
		printf("逆离散余弦变换输出文件打开失败\n");
		exit(-1);
	}
	unsigned char *idct_data=(unsigned char *)malloc(sizeof(char)*OffSet);
	fseek(fp_bmp,0L,SEEK_SET);
	fread(idct_data,sizeof(char),OffSet,fp_bmp);
	fwrite(idct_data,sizeof(char),OffSet,fp_idct_out);
	free(idct_data);
	unsigned char  *idct_pixel_out=(unsigned char *)malloc(sizeof(char)*line_data);
	for(row=0;row<bmpHeight;row++)
	{
		for(col=0;col<bmpWidth;col++)
		{
			idct_pixel_out[col*3+2]=(unsigned char)idct_r_array[bmpHeight-1-row][col];
			idct_pixel_out[col*3+1]=(unsigned char)idct_g_array[bmpHeight-1-row][col];
			idct_pixel_out[col*3+0]=(unsigned char)idct_b_array[bmpHeight-1-row][col];
		}
		fwrite(idct_pixel_out,sizeof(char),line_data,fp_idct_out);
	}	
	

	my_free();//全局变量－二维数组释放
	free(fp_gray_data);
	free(fp_data);
	free(temp_bmp_arry);

	end_time=clock();
	printf("耗时:%lf秒\n",(end_time-start_time)/CLOCKS_PER_SEC);
	return 0;
}

