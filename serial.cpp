
 
#include<bits/stdc++.h>

// #include<time.h>
using namespace std;
#define ll long long int
#define pb push_back
#define mp make_pair
#define mod 1000000007
#define maxn 10000001
#define endl "\n"
#define inf INT_MAX
#define vec vector <ll>
#define mps map <ll,ll> 
#define ppl vector<pair<ll,ll>> 
#define fi first
#define si second
#define pr pair<ll,ll> 
#define mvec map<ll,vector <ll>>
#define all(v) v.begin(),v.end()
#define infi LLONG_MAX
#define pi 3.141592653
#define ninf INT_MIN
#define ninfi LLONG_MIN
#define lbn lower_bound
#define ubn upper_bound
#define memset(a,b) memset(a,(b),sizeof(a))


vector<double> array_fill(ll i, ll n, ll v) {
    vector<double> a;
    for (; i < n; i++) {
        a.pb(v);
    }
    return a;
}
void printsol(vector<vector<double>>&A)
{
int i, j;
int row=A.size();
int col=A[0].size();

  for( i=0; i<row; i++ ) {
    for( j=0; j<col-1; j++ ) {
      if( A[i][j]!= 0 ) {
        cout << A[i][j] << "x" << j;
        if( j<A.size()-1) cout << " + ";
      }
      else
        cout << "      ";
    }
    cout << " = " << A[i][col-1] << endl;
  }
}


vector<double> gauss(vector<vector<double>>&A, vector<double> &x) {

    ll i, k, j;

    // Just make a single matrix
    for (i=0; i < A.size(); i++) { 
        A[i].pb(x[i]);
    }
    

    ll n = A.size();
    
    for (i=0; i < n; i++) { 
        // Search for maximum in this column
        double maxEl = abs(A[i][i]),
            maxRow = i;
        for (k=i+1; k < n; k++) { 
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }

       
         

        


        // Swap maximum row with current row (column by column)
        for (k=i; k < n+1; k++) { 
            double tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

      



        

        // Make all rows below this one 0 in current column
        for (k=i+1; k < n; k++) { 

            double c = (-1)*A[k][i]/A[i][i];
            for (j=i; j < n+1; j++) { 
                if (i==j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }

      
        }

    

     

        

    }

  // for(int i=0;i<A.size();i++)
  //   {
  //     for(int j=0;j<A[0].size();j++)
  //       cout<<A[i][j]<<" ";
  //     cout<<endl;
  //   }

  printsol(A);




//     // Solve equation Ax=b for an upper triangular matrix A
    x = array_fill(0, n, 0);
    
    for (i=n-1; i > -1; i--) { 
        x[i] = A[i][n]/A[i][i];
        for (k=i-1; k > -1; k--) { 
            A[k][n] -= A[k][i] * x[i];
        }
    }

   

    return x;
}


void read(vector<vector<double>>&A,vector<double>&x,int n)
{
  
      
    double r;

    
    for(ll i=0; i<n; i++ ) {      // creates a matrix of random
    x[i] = 0.0;             // integers, assuring non-
    for( ll j=0; j<n; j++ ) {      // scalar equations
      r = rand();
      A[i][j] = r;
      x[i] += j*r;
    }
  }

    // for(ll i=0;i<n;i++)
    // {
    //     for(ll j=0;j<n;j++)
    //         cin>>A[i][j];
    // }

    // for(ll i=0;i<n;i++)
    //     cin>>x[i];

}



 
int main()
{
    time_t start, end; 
     
  
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    ll n;
    cin>>n;
    vector<double> x(n);
    vector<vector<double>> A(n,vector<double>(n));
    read(A,x,n);




  

   



    
    
   clock_t time_req;
   time_req = clock();


    
    vector<double> ans=gauss(A,x);

    time_req = clock() - time_req;
    
    
    cout<<"Total Time Taken:"<<(float)time_req/CLOCKS_PER_SEC<<endl;

    
    for(ll i=0;i<x.size();i++)
    	cout<<"x"<<"["<<i<<"]"<<"="<<x[i]<<endl;
    


    
    return 0;
}