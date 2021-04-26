#include <iostream>
#include<armadillo>

using namespace std;
using namespace arma;
vector<mat> triang_sup(mat a,colvec b)
{
    int n=a.n_rows;
    mat a_aumt=join_rows(a,b);
    for(int k=0;k<n-1;k++)
    {
        for(int i=k+1;i<n;i++)
        {
            double mik=a_aumt(i,k)/a_aumt(k,k);
            for(int j=k;j<n+1;j++)
            {
                a_aumt(i,j)=a_aumt(i,j)-mik*a_aumt(k,j);
            }
        }
    }
    vector<mat>out;
    out.push_back(a_aumt.cols(0,n-1));
    out.push_back(a_aumt.cols(n,n));
    return out;
}
mat sust_atras(mat a,mat b)
{
    int n=a.n_rows;
    mat xk=zeros(1,n);
    for(int i=n-1;i>-1;i--)
    {
        double suma=0;
        for(int j=i+1;j<n;j++)
        {
            suma+=a(i,j)*xk(0,j);
        }
        double xi=(1/a(i,i))*(b(i,0)-suma);
        xk(0,i)=xi;
    }
    return xk;
}
mat gaussiana(mat a,mat b)
{
    if(det(a)!=0)
    {
        vector<mat> t=triang_sup(a,b);
        return sust_atras(t[0],t[1]);
    }
    cout<<"No es posible aplicar el metodo, la matriz no tiene inversa"<<endl;
    return {{0,0}};
}
int main()
{
    mat A={{2,-6,12,16},
           {1,-2,6,6},
           {-1,3,-3,-7},
           {0,4,3,-6}};

    colvec B={70,26,-30,-26};
    mat xk=gaussiana(A,B);
    cout<<"Las soluciones encontradas son:\nxk= "<<xk<<endl;
    return 0;
}


