#include<bits/stdc++.h>
using namespace std;
#include <stdlib.h>
#define ll long long int

#include<mpi.h>

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

    if( A[k*n+k] != 0 )     // calculates new value
      y[k] = b[k] / A[k*n+k];   // for equation solution
    else
      y[k] = 0.0;

    A[k*n+k] = 1.0;     // sets UTM diagonal value

    for( i=k+1; i<n; i++ ) {    // Guassian elimination occurs
      for( j=k+1; j<n; j++ )    // in all remaining rows
        A[i*n+j] -= A[i*n+k] * A[k*n+j];

      b[i] -= A[i*n+k] * y[k];
      A[i*n+k] = 0.0;
    }
  }
  
    
  
  return tstart;    // returns timing of serial algorithm
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
  double tstart, tfinish, TotalTime;    // timing decls

  if( argc < 2 ) {
    cout << "Usage\n";
    cout << "  Arg1 = number of equations / unkowns\n";
    return -1;
  }



  n = atoi(argv[1]);

  A = new double[n*n];        // space for matricies
  b = new double[n];
  y = new double[n];

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
  MPI_Comm com = MPI_COMM_WORLD;

  int size,rank;                          // Get rank/size info
  MPI_Comm_size(com,&size);
  MPI_Comm_rank(com,&rank);

  int manager = (rank == 0);    
          // Set the manager flag
                                          // if on processor 0

                                           // go serial if only 1 PE
  if (size == 1)
  {
    
    tstart = serial_gaussian ( A, b, y, n);

    tfinish = MPI_Wtime();
    cout<<"TotalTime is"<<tfinish-tstart<<endl;


   }
  
  else
  {

  if ( ( n % size ) != 0 )
  {
    // cout<<size<<" "<<n<<endl;
  cout << "Unknowns must be multiple of processors." << endl;
  return -1;
  }

  int np = (int) n/size;
  a = new double[n*np];     // size of each submatrix
  tmp = new double[n*np];   // and temp. submatrix

  if ( manager )
  {
  tstart = MPI_Wtime();   // first PE starts timing
  final_y = new double[n];  // and collection var

  }
          
  MPI_Scatter(A,n*np,MPI_DOUBLE,a,n*np,MPI_DOUBLE,0,com);

      
  for ( i=0; i < (rank*np); i++ )
  {
  // Get data (blocking)
  MPI_Bcast(tmp,n,MPI_DOUBLE,i/np,com);      //Each Process is Broadcasting np times
  MPI_Bcast(&(y[i]),1,MPI_DOUBLE,i/np,com);

      
  for (row=0; row<np; row++)
  {
  for ( j=i+1; j<n; j++ )
    a[row*n+j] = a[row*n+j] - a[row*n+i]*tmp[j];
  b[rank*np+row] = b[rank*np+row] - a[row*n+i]*y[i];
  a[row*n+i] = 0;
     }
  }

        // after receiving data from all prior
        // rows, begin calculcations
for (row=0; row<np; row++)
{
  for ( j=rank*np+row+1; j < n ; j++ )
{
  a[row*n+j] = a[row*n+j] / a[row*n+np*rank+row];
}
  y[rank*np+row] = b[rank*np+row] / a[row*n+rank*np+row];
  a[row*n+rank*np+row] = 1;

        // send your row's calculated data
  for ( i=0; i<n ; i++ )
  tmp[i] = a[row*n+i];

  MPI_Bcast (tmp,n,MPI_DOUBLE,rank,com);
  MPI_Bcast (&(y[rank*np+row]),1,MPI_DOUBLE,rank,com);

  // update lower rows
  for ( i=row+1; i<np; i++)
  {
  for ( j=rank*np+row+1; j<n; j++ )
    a[i*n+j] = a[i*n+j] - a[i*n+row+rank*np]*tmp[j];
  b[rank*np+i] = b[rank*np+i] - a[i*n+row+rank*np]*y[rank*np+row];
  a[i*n+row+rank*np] = 0;
  }
}
        // set idle and receive lower row's
        // broadcast message. VERY BAD, this
        // is why we implement cyclic mapping
        // in the other version
  for (i=(rank+1)*np ; i<n ; i++)
  {
  MPI_Bcast (tmp,n,MPI_DOUBLE,i/np,com);
  MPI_Bcast (&(y[i]),1,MPI_DOUBLE,i/np,com);
  }

        // Synchronize and gather rows
        // and y values
  MPI_Barrier(com);
  MPI_Gather(a,n*np,MPI_DOUBLE,A,n*np,MPI_DOUBLE,0,com);
  MPI_Gather(&(y[rank*np]),np,MPI_DOUBLE,final_y,np,MPI_DOUBLE,0,com);

  if (manager || (size==1) )
  {
    
    tfinish = MPI_Wtime();  // set and output timing val
    TotalTime = tfinish - tstart;
    cout << "TotalTime is:"<<TotalTime;
    cout << endl;
    
  

    // print_equations( A, y, n );  
    // print_solution(A,y,n);

  }
  
  

  MPI_Finalize();

    
}         
  
}
