#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <time.h>
#include <omp.h>
using namespace std;
#define PI 3.14159265
double sinc(double x);
int main()
{
	//time_t start_time,end_time;
	double start,finish;
	start=omp_get_wtime();
	ifstream F("dataF.txt"),q("dataq.txt"),X("dataX.txt");
	double TestSin[60],ArrayF[214][60],Arrayq[1][60],ArrayX[214][3],outputQ[60],dist[214][214];
	int i,j,k,l;
	//start_time = clock();
	//for(i=0;i<214;i++)
	for(j=0;j<60;j++)
		for(i=0;i<214;i++)
	{
		
		F>>ArrayF[i][j];
	}

	for(j=0;j<60;j++)
	{
		q>>Arrayq[1][j];
	}
	for(j=0;j<3;j++)
		for(i=0;i<214;i++)
	{
		X>>ArrayX[i][j];
	}


	//#pragma omp parallel for private(j,l)
	for(i=0;i<214;i++)
		for(j=0;j<214;j++)
		{
			dist[i][j]=0;
			for(l=0;l<3;l++)
			{	
				dist[i][j]+=(ArrayX[i][l]-ArrayX[j][l])*(ArrayX[i][l]-ArrayX[j][l]);
			}
			dist[i][j]=sqrt(dist[i][j]);
			
		}

	//cout<<sin(30)<<endl;
	double sum;

#pragma omp parallel num_threads(60) private(i,j,k) 
{
	
	k=omp_get_thread_num();
	{
		outputQ[k]=0;
		
		for(i=0;i<214;i++)
		{

			for(j=0;j<214;j++)
			{			
				outputQ[k] = outputQ[k] + ArrayF[i][k]*ArrayF[j][k]*sinc(Arrayq[1][k]*dist[i][j]);
			}
		}

	}
	

}
for(i=0;i<60;i++)cout<<outputQ[i]<<endl;
	//end_time = clock();
	//cout<<" total time is "<<(double)(end_time-start_time)/CLOCKS_PER_SEC<<" sec"<<endl;
	finish = omp_get_wtime();
	cout<<"total time is "<<(finish-start)<<endl;
}	

double sinc(double x)
{
	if(x==0)
		return 1;
	else return sin(x)/x;
}	
