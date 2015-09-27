//Distributed two-dimensionnal fft
//Xiaofei Qiu
//ECE 8893 Project 1



#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <signal.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

#include "Complex.h"
#include "InputImage.h"
using namespace std;

int nCpus, rank, rc;
Complex* Garbage;


void Transform1D(Complex* h, int N, Complex* H);
void Transform2D(const char* inputFN);
void CheckReturnCode(const int& rc);
void Transpose(Complex *source,const int& rows, const int& cols);


int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  rc = MPI_Init(&argc,&argv);
  Transform2D(fn.c_str()); // Perform the transform.

  delete [] Garbage;


   MPI_Finalize();
}

////////////////////////////////////////////////
//                 Routines                  //

void Transform2D(const char* inputFN)
{
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// MPI setup
    MPI_Comm_size(MPI_COMM_WORLD,&nCpus);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Request request;
    MPI_Status status;

// Get source data
    InputImage image(inputFN);
    Complex *Data = image.GetImageData();
    Garbage = Data;
// Const setup
    const int nCols = image.GetWidth();
    const int nRows = image.GetHeight();
    const int rRound = nRows / nCpus;
    const int cRound = nCols / nCpus;
    const int startRow = rRound * rank;
//    double masterBuffer[nRows*nCols*2];

// Output container setup
    Complex rOut[nCols*rRound];
    Complex cOut[nRows*cRound];
    Complex cIn[nRows*cRound];


/////////////////////////////////////////////////////////////////////////
//                DFT on rows
/////////////////////////////////////////////////////////////////////////

    for(int i=0;i<rRound;i++)
    {
        Transform1D(&Data[(startRow+i)*nCols],nCols,&rOut[i*nCols]);
    }

//////////////////////////////////////////////////////////////////////////
//                 Send and recv data
//////////////////////////////////////////////////////////////////////////

    if(rank != 0)
    {
        rc = MPI_Send(rOut,sizeof(rOut),MPI_CHAR,0,0,MPI_COMM_WORLD);
        CheckReturnCode(rc);
    }

    if(rank == 0)
    {
        for(int j=1; j<nCpus;j++)
        {
            rc = MPI_Recv(&Data[rRound*j*nCols],sizeof(rOut),MPI_CHAR,j,0,MPI_COMM_WORLD,&status);
            CheckReturnCode(rc);
        }

        for(int i=0;i<nRows*cRound;i++)
        {
            Data[i] = rOut[i];
        }

        Transpose(Data,nRows,nCols);

        for(int i=0;i<nCpus;i++)
        {
            rc = MPI_Isend(&Data[i*cRound*nRows],sizeof(cIn),MPI_CHAR,i,0,MPI_COMM_WORLD,&request);
            CheckReturnCode(rc);
        }
    }

    MPI_Irecv(cIn,sizeof(cIn),MPI_CHAR,0,0,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,&status);

//////////////////////////////////////////////////////////////////////////
//                 DFT on cols
/////////////////////////////////////////////////////////////////////////

    for(int i=0;i<cRound;i++)
    {
        Transform1D(&cIn[i*nRows],nRows,&cOut[i*nRows]);
    }

/////////////////////////////////////////////////////////////////////////
//                 Send to Master
/////////////////////////////////////////////////////////////////////////

    if(rank!=0)
    {
        MPI_Send(cOut,sizeof(cOut),MPI_CHAR,0,0,MPI_COMM_WORLD);
    }

    if(rank==0)
    {
        for(int i=1;i<nCpus;i++)
        {
            rc = MPI_Recv(&Data[i*cRound*nCols],sizeof(cOut),MPI_CHAR,i,0,MPI_COMM_WORLD,&status);
            CheckReturnCode(rc);
        }

        for(int i=0;i<nRows*cRound;i++)
        {
            Data[i] = cOut[i];
        }

        Transpose(Data,nRows,nCols);
        image.SaveImageData("MyAfter2d.txt",Data,nCols,nRows);
    }
}

void Transpose(Complex* source,const int& rows, const int& cols)
{
    Complex *temp = new Complex[cols*rows];
    int r = 0; int c = 0; int index = 0;

    for(int i =0;i<cols*rows;i++)
    {
        temp[i] = source[i];
    }

    for(int j = 0; j<cols*rows;j++)
    {
        r = j/cols;
        c = j%cols;
        index = c*rows + r;
        source[index] = temp[j];
    }
    delete [] temp;
}

void Transform1D(Complex* h, int N, Complex* H)
{
    Complex W;
    Complex sum(0,0);
    for(int n=0;n<N;n++)
    {
        for(int k = 0;k<N;k++)
        {
            W.real = cos(2*M_PI*n*k/N);
            W.imag = -sin(2*M_PI*n*k/N);
            sum = sum + W * h[k];
        }
        H[n] = sum;
        sum.real=0;
        sum.imag=0;
    }
}

void CheckReturnCode(const int& rc)
{
        if (rc != MPI_SUCCESS)
            {
              cout << "Rank " << rank
                   << " send failed, rc " << rc << endl;
              MPI_Finalize();
              exit(1);
            }
}
