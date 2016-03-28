#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <time.h>
#include <mpi.h>
//need 11 processors
using namespace std;
#define PI 3.14159265
double sinc(double x);
int main(int argc,char *argv[])
{
    char processor_name[MPI_MAX_PROCESSOR_NAME];	
    double start_time,end_time;
    int numprocessor,rank,namelen;
    int tag = 99;
    ifstream F("dataF.txt"),q("dataq.txt"),X("dataX.txt");
    double TestSin[60],ArrayF[214][60],Arrayq[60],ArrayX[214][3],outputQ[60],dist[214][214];
    int i,j,k,l;
    int step;
    //int output_start=0;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocessor);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Get_processor_name(processor_name,&namelen);    
//	ifstream F("dataF.txt"),q("dataq.txt"),X("dataX.txt");
//	double TestSin[60],ArrayF[214][60],Arrayq[1][60],ArrayX[214][3],outputQ[60],dist[214][214];
//	int i,j,k,l;
    MPI_Barrier(MPI_COMM_WORLD);
	start_time = MPI_Wtime();
    step=60/(numprocessor-1);
    //parallel read data
    MPI_Status status;
    //read X used for calculate
    if(rank==0)
    {
        for(j=0;j<3;j++)
            for(i=0;i<214;i++)
            {
                X>>ArrayX[i][j];
            }
        for(j=0;j<60;j++)
            for(i=0;i<214;i++)
            {
                
                F>>ArrayF[i][j];
            }
        
        for(j=0;j<60;j++)
        {
            q>>Arrayq[j];
	    //cout<<Arrayq[1][j]<<endl;
        }
        
    }
	MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(ArrayX,3*214,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(Arrayq,60,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(ArrayF,60*214,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	//use 20 processor to calculate
	if(rank==1)
{
/*
for(j=0;j<3;j++)
   {         for(i=0;i<214;i++)
            {
                cout<<ArrayX[i][j]<<" ";
            }
cout<<endl;
}


        for(j=0;j<2;j++)
          {  for(i=0;i<214;i++)
            {

                cout<<ArrayF[i][j]<<" ";;
            }
cout<<endl;
}

 for(j=0;j<60;j++)
        {
            //q>>Arrayq[1][j];
            cout<<"ll"<<Arrayq[j]<<endl;
                    }
*/
}
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
    {

        for(i=1;i<numprocessor;i++)
        {
            
            MPI_Recv(&outputQ[(i-1)*step],step,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
            for(j=0;j<step;j++)
                cout<< outputQ[j+(i-1)*step]<<endl;
        }
   

     }
    else
    {
	
        int output_start=(rank-1)*step;
        for(k=output_start;k<output_start+step;k++)
        {
            outputQ[k]=0;
            for(i=0;i<214;i++)
            {
                for(j=0;j<214;j++)
                {
                    //TestSin[k]+=Arrayq[1][k]*dist[i][j];
                    dist[i][j]=0;
                    for(l=0;l<3;l++)
                    {
                        dist[i][j]+=(ArrayX[i][l]-ArrayX[j][l])*(ArrayX[i][l]-ArrayX[j][l]);
                    }
                    dist[i][j]=sqrt(dist[i][j]);
                    outputQ[k] = outputQ[k] + ArrayF[i][k]*ArrayF[j][k]*sinc(Arrayq[k]*dist[i][j]);
                }
            }
       //    cout<< outputQ[k]<<endl;
           //printf("output %f from %d\n",outputQ[k],rank);
        }
        MPI_Send(&outputQ[output_start],step,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        
    }

    MPI_Barrier(MPI_COMM_WORLD);

	end_time = MPI_Wtime();
    //MPI_Finalize();
    if(rank==0)
    {
        cout<<" total time is "<<(end_time-start_time)<<" sec"<<endl;
        
    }
	
    MPI_Finalize();

}

double sinc(double x)
{
	if(x==0)
		return 1;
	else return sin(x)/x;
}	
