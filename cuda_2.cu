#include <iostream>
#include <fstream>
#include <cuda.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

#define PI 3.14159265
__host__ void readdata(double **ArrayX,double **ArrayF,double *Arrayq)
{
    
    
}
__device__ double sinc(double x)
{
	if(x==0)
		return 1;
	else return sin(x)/x;
}
//since max number of threads is 1024 per block which is smaller than 214*214
//so separate the calculation,214 x blocks and 214 threads for each of 60 y blocks
//by which I can use shared memory to do reducation
 
__global__ void cal_output_first(double *ArrayF,double *Arrayq,double *dist,double *output_first)
{
    int i,j,k;
    i=threadIdx.x;
    j=blockIdx.x;
    k=blockIdx.y;
    __shared__ double sdata[214];

    int tid=threadIdx.x;
	if(threadIdx.x<214&&blockIdx.y<60&&blockIdx.x<214)
{
	sdata[tid]=ArrayF[i+214*k]*ArrayF[j+214*k]*sinc(Arrayq[k]*dist[i+214*j]);
    __syncthreads();
    double sum=0;

    for(int m=0;m<214;m++)
    {
        sum+=sdata[m];
        __syncthreads();
    }
    if(tid==0)
	output_first[blockIdx.x+blockIdx.y*blockDim.x]=sum;
}
}
//calculate sum of remaining 214 value for outputq[60]
//214 threads and 60 x blocks
__global__ void cal_output_second(double *output_first,double *outputQ)
{
    int i,k;
    i=threadIdx.x;
    k=blockIdx.x;
    __shared__ double sdata[214]; 
if(threadIdx.x<214&&blockIdx.x<60)
{
    int tid=threadIdx.x;
        sdata[tid]=output_first[i+214*k];
    __syncthreads();
    double sum=0;

    for(int m=0;m<214;m++)
    {
        sum+=sdata[m];
        __syncthreads();
    }
    if(tid==0)
        outputQ[blockIdx.x]=sum;
}
}
__global__ void cal_dist(double *ArrayX,double *dist)
{
    int i,j;
    i=threadIdx.x;
    j=blockIdx.x;
if(i<214&&j<214)
{
    double temp1,temp2,temp3;
    temp1=(ArrayX[i]-ArrayX[j])*(ArrayX[i]-ArrayX[j]);
    temp2=(ArrayX[i+214]-ArrayX[j+214])*(ArrayX[i+214]-ArrayX[j+214]);
    temp3=(ArrayX[i+428]-ArrayX[j+428])*(ArrayX[i+428]-ArrayX[j+428]);
	
    dist[i+214*j]=sqrt(temp1+temp2+temp3);

//printf("temp is %d",temp);
}
}
int main()
{
    cudaEvent_t start = 0;
    cudaEvent_t stop = 0;
    float time = 0;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    double ArrayF[214*60],Arrayq[60],ArrayX[214*3],outputQ[60];
    int i,j;
    int sizeF=214*60,sizeq=60,sizeX=214*3;
    int sizeQ=60,sizedist=214*214;
    ifstream F("dataF.txt"),q("dataq.txt"),X("dataX.txt");
   
    for(j=0;j<3;j++)
        for(i=0;i<214;i++)
        {
            X>>ArrayX[i+j*214];
        }
    for(j=0;j<60;j++)
        for(i=0;i<214;i++)
        {
            
            F>>ArrayF[i+214*j];
        }
    
    for(j=0;j<60;j++)
    {
        q>>Arrayq[j];
        //cout<<Arrayq[1][j]<<endl;
    }
    //calculate distance dist[214*214]
    int dsizeX=sizeX*sizeof(double);
    int dsizedist=sizedist*sizeof(double);
    int dsizeF=sizeF*sizeof(double);
    int dsizeq=sizeq*sizeof(double);
    int dsizeQ=sizeQ*sizeof(double);
    int dsize_first=214*60*sizeof(double);
    double *d_dist,*d_ArrayX,*d_ArrayF,*d_Arrayq,*d_outputQ;
    double *d_output_first;
    //cudaEventRecord(start,0);
    cudaMalloc((void**)&d_dist,dsizedist);
    cudaMalloc((void**)&d_ArrayX,dsizeX);
    cudaMalloc((void**)&d_ArrayF,dsizeF);
    cudaMalloc((void**)&d_Arrayq,dsizeq);
    cudaMalloc((void**)&d_outputQ,dsizeQ);
    //allocate memory for output_first
    cudaMalloc((void**)&d_output_first,dsize_first);
    cudaMemcpy(d_ArrayX,&ArrayX,dsizeX,cudaMemcpyHostToDevice);
    cudaMemcpy(d_ArrayF,&ArrayF,dsizeF,cudaMemcpyHostToDevice);
    cudaMemcpy(d_Arrayq,&Arrayq,dsizeq,cudaMemcpyHostToDevice);
    cudaMemcpy(d_outputQ,&outputQ,dsizeQ,cudaMemcpyHostToDevice);
    dim3 DimGrid1(256,1,1);
    dim3 DimBlock1(256, 1,1);
    dim3 DimGrid2(256,64,1);
    dim3 DimBlock2(256, 1,1);
    dim3 DimGrid3(64,1,1);
    dim3 DimBlock3(256, 1,1);
    cudaEventRecord(start,0);
    cal_dist<<<DimGrid1,DimBlock1>>>(d_ArrayX,d_dist);
    //cudaEventRecord(stop,0);
    //cudaEventSynchronize(stop);
    cudaThreadSynchronize();
    cal_output_first<<<DimGrid2,DimBlock2>>>(d_ArrayF,d_Arrayq,d_dist,d_output_first);
    //output_first size 214*60
    //calculate output
    cal_output_second<<<DimGrid3,DimBlock3>>>(d_output_first,d_outputQ);
    //cudaMemcpy(&outputQ,d_outputQ,dsizeQ,cudaMemcpyDeviceToHost);
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaMemcpy(&outputQ,d_outputQ,dsizeQ,cudaMemcpyDeviceToHost);
    cudaFree(d_dist);
    cudaFree(d_ArrayX);
    cudaFree(d_ArrayF);
    cudaFree(d_Arrayq);
    cudaFree(d_outputQ);
    //Free memory for output_first
    cudaFree(d_output_first);

    for(i=0;i<60;i++)
    {
	cout<<outputQ[i]<<endl;
    }
    /*
    for(j=0;j<214;j++)
	{
		for(i=0;i<214;i++)
	cout<<i<<","<<j<<"__"<<dist[i+j*214]<<endl;
	cout<<endl;
	}
    */
    	cudaEventElapsedTime(&time,start,stop);
        cout<<"Time for the kernel: "<<time<<endl;
    return EXIT_SUCCESS;
     
}

