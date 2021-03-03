

#include <iostream>
#include<bits/stdc++.h>
#include <stdlib.h>

#include "mpi.h"

using namespace std;
// Serial Algorithm included for easy comparative timings

double serial_gaussian( double *A, double *b, double *y, int n )
{
  int i, j, k;
  double tstart = MPI_Wtime();

  for( k=0; k<n; k++ ) {    // k = current row
    for( j=k+1; j<n; j++ ) {    // in division step
      if( A[k*n+k] != 0)
        A[k*n+j] = A[k*n+j] / A[k*n+k];
      else
        A[k*n+j] = 0;
    }

    if( A[k*n+k] != 0 )     // calc new equation
      y[k] = b[k] / A[k*n+k];   // value
    else
      y[k] = 0.0;

    A[k*n+k] = 1.0;     // set diagonal of UTM

    for( i=k+1; i<n; i++ ) {    // update all preceeding rows
      for( j=k+1; j<n; j++ )
        A[i*n+j] -= A[i*n+k] * A[k*n+j];

      b[i] -= A[i*n+k] * y[k];
      A[i*n+k] = 0.0;
    }
  }
  return tstart;      // returns timing value
}



void print_equations( double *A, double *y, int n )
{
  int i, j;

  for( i=0; i<n; i++ ) {
    for( j=0; j<n; j++ ) {
      if( A[i*n+j] != 0 ) {
        cout << A[i*n+j] << "x" << j;
        if( j<n-1 ) cout << " + ";
      }
      else
        cout << "      ";
    }
    cout << " = " << y[i] << endl;
  }
}
void print_solution( double *A, double *y, int n )
{
  int i, j;
  double *x = new double[n];

  for( i=n-1; i>=0; i-- ) {
    x[i] = y[i];
    for( j=n-1; j>i; j-- )
      x[i] -= x[j] * A[i*n+j];
  }
  for( i=0; i<n; i+=4 ) {
    for( j=i; j<i+4 && j<n; j++ )
      cout << "x[" << j << "] = " << x[j] << "  ";
    cout << endl;
  }
}
int main( int argc, char *argv[] )
{
  double *A, *b, *y, *a, *tmp, *final_y;  // var decls
  int i, j, n, row, r;
  double tstart,tfinish,TotalTime;    // timing vars

  if( argc < 2 ) {
    cout << "Usage\n";
    cout << "  Arg1 = number of equations / unkowns\n";
    return -1;
  }

  n = atoi(argv[1]);

  A = new double[n*n];        // space for matricies
  b = new double[n];
  y = new double[n];
  tmp = new double[n];

/* Coefficient computer: */
  for( i=0; i<n; i++ ) {    // creates a matrix of random
    b[i] = 0.0;       // integers, assuring non-
    for( j=0; j<n; j++ ) {    // scalar equations
      r = rand();
      A[i*n+j] = r;
      b[i] += j*r;
    }
  }

  MPI_Init (&argc,&argv);   // Initialize MPI
  MPI_Status status;
  MPI_Comm com = MPI_COMM_WORLD;

  int size,rank;                          // Get rank/size info
  MPI_Comm_size(com,&size);
  MPI_Comm_rank(com,&rank);

  int manager = (rank == 0);              // Set the manager flag
                                          // if on processor 0

            // go serial if only 1 PE
  if (size == 1)
  tstart = serial_gaussian ( A, b, y, n);
  else
  {

          // striped mapping
          // allows only p < i*n
    if ( ( n % size ) != 0 )
    {
      cout << "Unknowns must be multiple of processors." << endl;
      return -1;
    }

    int np = (int) n/size;
    a = new double[n*np];     // size of each submatrix

    if ( manager )
    {
    tstart = MPI_Wtime();   // first PE starts timing
    final_y = new double[n];  // and collection var

    }

          // divys up the rows to PEs
          // in blocked mapping fashion
    MPI_Scatter(A,n*np,MPI_DOUBLE,a,n*np,MPI_DOUBLE,0,com);

       
    for ( i=0; i < rank*np; i++ )
    {
          // Get data (blocking)
      MPI_Recv(tmp,n,MPI_DOUBLE,rank-1,20,com,&status);
      MPI_Recv(&(y[i]),1,MPI_DOUBLE,rank-1,30,com,&status);

          
      if ( rank != (size-1) )
      {
        MPI_Send (tmp,n,MPI_DOUBLE,rank+1,20,com);
        MPI_Send (&(y[i]),1,MPI_DOUBLE,rank+1,30,com);
      }

       
    for (row=0; row<np; row++)
    {
      for ( j=i+1; j<n; j++ )
        a[row*n+j] = a[row*n+j] - a[row*n+i]*tmp[j];
        b[rank*np+row] = b[rank*np+row] - a[row*n+i]*y[i];
        a[row*n+i] = 0;
      }
    }

      
   for (row=0; row<np; row++)
   {
      for ( j=rank*np+row+1; j < n ; j++ )
        a[row*n+j] = a[row*n+j] / a[row*n+np*rank+row];
      y[rank*np+row] = b[rank*np+row] / a[row*n+np*rank+row];
      a[row*n+np*rank+row] = 1;

      for (i=0; i<n; i++)
        tmp[i] = a[row*n+i];

          
        if ( rank != (size-1) )
        {
          MPI_Send (tmp,n,MPI_DOUBLE,rank+1,20,com);
          MPI_Send (&(y[rank*np+row]),1,MPI_DOUBLE,rank+1,30,com);
        }

  // Update lower rows

    for (i=row+1; i<np; i++)
    {
      for (j=rank*np+row+1; j<n; j++)
        a[i*n+j] = a[i*n+j] - a[i*n+rank*np+row]*tmp[j];
      b[rank*np+i] = b[rank*np+i] - a[i*n+rank*np+row]*y[rank*np+row];
      a[i*n+rank*np+row] = 0;
    }
  }
        // Synchronize and gather rows
        // and y values
  MPI_Barrier(com);
  MPI_Gather(a,n*np,MPI_DOUBLE,A,n*np,MPI_DOUBLE,0,com);
  MPI_Gather(&(y[rank*np]),np,MPI_DOUBLE,final_y,np,MPI_DOUBLE,0,com);
  //MPI_Finalize();

  y = final_y;      // reconstruct y vector with gathered
        // values

  }         // end of ELSE
  if (manager || (size==1) )
  {
    tfinish = MPI_Wtime();
    TotalTime = tfinish - tstart;

    // print_equations(A,y,n);
    // print_solution(A,y,n);
    cout << "TIME TAKEN  : " << TotalTime;    // Output Time
    cout << endl;
//    print_equations( A, y, n ); // Output matrix (not wanted in timings)
  }

  MPI_Finalize();
}
