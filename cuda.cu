#include <iostream>
#include <cuda.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

#define PI 3.14159265
double sinc(double x);
__host__ void inital(float *v1,float *v2,int vsize)
{
	srand(time(NULL));
	int i;
	for(i=0;i<vsize;i++)
	{
		v1[i]=static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/10));	

		v2[i]=static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/10));
		//cout<<v1[i]<<"--"<<v2[i]<<endl;
	}
}
__global__ void VecAdd(float *v1,float *v2,float *v3,int n)
{
	int i = threadIdx.x+blockDim.x*blockIdx.x+blockDim.x*blockDim.x*blockIdx.y;
	//cout<<i<<endl;
	if(i<n)
	{
		v3[i]=v1[i]+v2[i];
	}
}

/*
int VecAdd(float *v1,float *v2,float *v3,int n)
{
	// Run ceil(n/1000) blocks of 1000 threads each
	dim3 DimGrid(ceil(n/1000.0), 1, 1);
	dim3 DimBlock(1000, 1, 1);
	VecAdd<<<DimGrid,DimBlock>>>(v1, v2, v3);	
}


__device__ float compute(float *v1,float *v2)
{

}
*/
int main()
{
	cudaEvent_t start = 0;
	cudaEvent_t stop = 0;
	float time = 0;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	double TestSin[60],ArrayF[214][60],Arrayq[60],ArrayX[214][3],outputQ[60],dist[214][214];

	//initialize vector
	inital(v1,v2,vsize);
	//check data
	/*
	for(i=0;i<10;i++)
	{
		cout<<v1[i]<<"--"<<v2[i]<<endl;
	}
	*/
 	//allocate memory 1000 threads per block
	int fsize = vsize*sizeof(float);
	
	float *v1_d,*v2_d,*v3_d;
	cudaMalloc((void**)&v1_d,fsize);
	cudaMemcpy(v1_d,v1,fsize,cudaMemcpyHostToDevice);
	cudaMalloc((void**)&v2_d,fsize);
        cudaMemcpy(v2_d,v2,fsize,cudaMemcpyHostToDevice);
	cudaMalloc((void**)&v3_d,fsize);
        
	//cudaMemcpy(v3_d,v3,fsize,cudaMemcpyHostToDevice);

	//kernel code
	//int bsize=vsize/1000;
	int gridsizex=ceil(vsize/1024.0);
	if(ceil(vsize/1024.0)>1024)
		gridsizex = 1024;
	
	int gridsizey = ceil(vsize/1024.0/1024.0);
	//cout<<gridsizex<<endl;
	dim3 DimGrid(gridsizex,gridsizey,1);
	dim3 DimBlock(1024, 1,1);
	cudaEventRecord(start,0);
	VecAdd<<<DimGrid,DimBlock>>>(v1_d,v2_d,v3_d,vsize);
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	//Synchronize
	cudaThreadSynchronize();
	//Transfer v3 back to host
	cudaMemcpy(v3,v3_d,fsize,cudaMemcpyDeviceToHost);
	
	//output
	//check	
	for(i=0;i<vsize;i++)
	{
		if(v3[i]!=v1[i]+v2[i])
		{
			cout<<i<<"--"<<v1[i]<<"--"<<v2[i]<<"--"<<v3[i]<<endl;
			break;
		}
	}
	
	//Free memory
	cudaFree(v1_d);cudaFree(v2_d);cudaFree(v3_d);

	cudaEventElapsedTime(&time,start,stop);
	cout<<"Time for the kernel: "<<time<<endl;
	return EXIT_SUCCESS;
}
