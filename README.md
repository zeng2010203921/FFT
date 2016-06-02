# FFT
测试图片：lena.bmp 512*512 真彩色
图片要求：长＝宽＝pow(2,n)　真彩色图片
解决办法：灰度图像，存在颜色索引表，像素值是颜色索引值
          32位图像，存在alpha通道
长！＝宽：像素值的填充

deal_image.c:实现图像的灰度变换、拉普拉斯变换，离散余弦变化
fft.c：一维实现图像的快速傅里叶变化
fft_dimen.c:图像的二维傅里叶变换
thread_fft.c:使用三线程单独处理RGB三方数据
cpu_affinity.c:线程亲和性，在某一CPU上运行
six_step_fft.c:六步快速傅里叶算法

测试工具：多核模拟器－sniper simulator
