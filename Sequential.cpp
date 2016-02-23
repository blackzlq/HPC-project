#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <time.h>

using namespace std;
#define PI 3.14159265
double sinc(double x);
int main()
{
	time_t start_time,end_time;
	ifstream F("dataF.txt"),q("dataq.txt"),X("dataX.txt");
	double TestSin[60],ArrayF[214][60],Arrayq[1][60],ArrayX[214][3],outputQ[60],dist[214][214];
	int i,j,k,l;
	start_time = clock();
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
	for(k=0;k<60;k++)
	{
		outputQ[k]=0;
		for(i=0;i<214;i++)
		{
			for(j=0;j<214;j++)
			{	
				//TestSin[k]+=Arrayq[1][k]*dist[i][j];
				outputQ[k] = outputQ[k] + ArrayF[i][k]*ArrayF[j][k]*sinc(Arrayq[1][k]*dist[i][j]);
			}
		}
		cout<< outputQ[k]<<endl;
	}
	end_time = clock();
	cout<<" total time is "<<(double)(end_time-start_time)/CLOCKS_PER_SEC<<" sec"<<endl;

}	

double sinc(double x)
{
	if(x==0)
		return 1;
	else return sin(x)/x;
}	
