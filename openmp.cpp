#include<bits/stdc++.h>
#include<omp.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
using namespace std;
#define ll long long int

#define SIZE 1200
#define MAXTHREADS 20
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





int main()
{

  double *A, *b, *y;
  A = new double[SIZE*SIZE];  
  int i,j;    
  b = new double[SIZE];
  y = new double[SIZE];
double r;

    
    for(ll i=0; i<SIZE; i++ ) {      // creates a matrix of random
    b[i] = 0.0;             // integers, assuring non-
    for( ll j=0; j<SIZE; j++ ) {      // scalar equations
      r = rand();
      A[i*SIZE+j] = r;
      b[i] += j*r;
    }
  }
int chunk, start, end,k,tid;
double start_time, end_time;
double diff_time;



cout<<"Size is "<<SIZE<<endl;


cout<<"AFTER THE READ FUNCTION"<<endl;
start_time =omp_get_wtime();
chunk = SIZE/MAXTHREADS;
omp_set_num_threads (MAXTHREADS) ;

#pragma omp parallel shared(A,b,y,i,j,chunk) private(k,tid,start,end)
{
  tid = omp_get_thread_num();
   start = tid*chunk;
end = start + chunk;
// printf("THREAD %d CREATED WITH START : %d AND END:%d \n",tid, start, end);

for(k=0;k<SIZE;k++)
{
if((k>=start)&&(k<end))
{
for(j = k+1;j<SIZE; j++)
{

  if(A[k*SIZE+k] !=0)
    A[k*SIZE+j] = A[k*SIZE+j] / A[k*SIZE+k];
  else
    A[k*SIZE+j] = 0;
}

if( A[k*SIZE+k] != 0)     // calculates new value
           y[k] = b[k] / A[k*SIZE+k];   // for equation solution
        else
            y[k] = 0.0;

        A[k*SIZE+k] = 1.0;


}
// #pragma omp barrier
#pragma omp for schedule(dynamic, chunk) nowait
for(i = k+1;i<SIZE;i++)
{
for(j =k+1;j<SIZE; j++)
A[i*SIZE+j] -= A[i*SIZE+k] * A[k*SIZE+j];
b[i] -= A[i*SIZE+k] * y[k];
A[i*SIZE+k] = 0.0;

}
// #pragma omp barrier

}
}// end of omp parallel.
// print_equations( A, y, SIZE);  // output matrix (not wanted in timings)
// print_solution(A,y,SIZE);
end_time = omp_get_wtime();
diff_time =end_time-start_time;
 


printf("TOTAL TIME : %f",diff_time);
printf("\n");









return 0;
}






