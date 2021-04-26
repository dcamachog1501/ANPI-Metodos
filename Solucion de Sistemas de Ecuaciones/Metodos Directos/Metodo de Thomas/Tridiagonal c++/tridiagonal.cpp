#include <iostream>
#include<armadillo>

using namespace std;
using namespace arma;

mat tridiagonal(double a,double b,double c,int n)
{
    mat out=zeros(n,n);
    for(int i=0;i<n;i++)
    {
        if(i==0)
        {
            out(i,i+1)=b;
        }
        else if(i==n-1)
        {
            out(i,i-1)=c;
        }
        else
        {
            out(i,i+1)=b;
            out(i,i-1)=c;
        }
        out(i,i)=a;
    }
    return out;
}
int main()
{
   cout<<tridiagonal(1,2,3,4)<<endl;
    return 0;
}


