#include <iostream>
#include<armadillo>

using namespace std;
using namespace arma;

bool equals(mat a,mat b)
{
    for(int i=0;i<a.n_rows;i++)
    {
        for(int j=0;j<a.n_cols;j++)
        {
            if(a(i,j)!=b(i,j))
            {
                return false;
            }
        }
    }
    return true;
}
bool is_symetric(mat a)
{
    int n=a.n_rows;
    int m=a.n_cols;
    mat a_t=a.t();
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            if(a(i,j)!=a_t(i,j))
            {
                return false;
            }
        }
    }
    return true;
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
mat sust_adelante(mat a,mat b)
{
    int n=a.n_rows;
    mat xk=zeros(1,n);
    for(int i=0;i<n;i++)
    {
        double suma=0;
        for(int j=0;j<i;j++)
        {
            suma+=a(i,j)*xk(0,j);
        }
        double xi=(1/a(i,i))*(b(i,0)-suma);
        xk(0,i)=xi;
    }
    return xk;
}
mat cholesky_l(mat a)
{
    int n=a.n_rows;
    int m=a.n_cols;
    if(n==m)
    {
        for(int i=0;i<n;i++)
        {
            mat temp=a.submat(0,0,i,i);
            if(det(temp)<0)
            {
                cout<<"Uno de los determinantes de las submatrices principales de a no es positivo"<<endl;
                return a;
            }
        }
        if(!is_symetric(a))
        {
            cout<<"La matriz provista no es simetrica"<<endl;
            return a;
        }
        mat l=zeros(n,m);
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<i+1;j++)
            {
                double sum=0;
                if(j==i)
                {
                    for(int k=0;k<j;k++)
                    {
                        sum+= pow(l(j, k), 2);
                    }
                    l(j,j)=sqrt(a(j,j)-sum);
                }
                else
                {
                    for(int k=0;k<j;k++)
                    {
                        sum+=l(i,k)*l(j,k);
                    }
                    l(i,j)=(a(i,j)-sum)/l(j,j);
                }
            }
        }
        return l;
    }
    cout<<"La matriz no es cuadrada"<<endl;
    return a;
}
mat fact_cholesky(mat a,mat b)
{
    mat l=cholesky_l(a);
    if(!equals(a,l))
    {
        mat yk=sust_adelante(l,b);
        return sust_atras(l.t(),yk.t());
    }
    cout<<"Aplicando el metodo de substitucion por A techo...\n"<<endl;
    mat at=a.t()*a;
    mat bt=a.t()*b;
    return fact_cholesky(at,bt);


}
int main()
{
    mat A={{1,-1,2},
           {-2,0,4},
           {0,-2,7}};

    colvec B={0,2,5};
    mat xk=fact_cholesky(A,B);
    cout<<"Las soluciones encontradas son:\nxk= "<<xk<<endl;
    return 0;
}


