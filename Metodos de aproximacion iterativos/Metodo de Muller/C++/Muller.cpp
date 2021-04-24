#include <iostream>
#include<ginac/ginac.h>
#include "cmath"

using namespace std;
using namespace GiNaC;

ex funct(string f,ex value)
{
    symbol x;
    symtab table;
    table["x"]=x;
    parser reader(table);
    ex expresion=reader(f);
    return evalf(expresion.subs(x==value));
}
vector<ex>muller(string f,ex x0,ex x1,ex x2,ex tol,int iterMax)
{
    vector<ex> output;
    ex xk = x1;
    ex error = 0;
    int k = 0;
    while (k < iterMax)
    {
        if(x0-x1==0||x0-x2==0||x1-x2==0)
        {
            cout<<"Indefinicion en alguna de las iteraciones,se retornan los valores hasta dicha iteracion"<<endl;
            break;
        }
        ex a=((x1-x2)*(funct(f,x0)-funct(f,x2))-(x0-x2)*(funct(f,x1)-funct(f,x2)))/((x0-x1)*(x0-x2)*(x1-x2));
        ex b=((pow((x0-x2),2))*(funct(f,x1)-funct(f,x2))-(pow((x1-x2),2))*(funct(f,x0)-funct(f,x2)))/((x0-x1)*(x0-x2)*(x1-x2));
        ex c= funct(f,x2);

        xk=x2-((2*c)/(b+(b/abs(b))* sqrt(pow(b,2)-4*a*c)));
        error=funct(f,xk);

        if (abs(error)<= tol)
        {
            break;
        }
        if(abs(x0-xk)<abs(x2-xk))
        {
            x2=x1;
        }
        else
        {
            x0=x1;
        }
        x1=xk;
        k++;
    }
    output.push_back(xk);
    output.push_back(abs(error));
    output.push_back(k);
    return output;
}

int main()
{
    Digits=10;
    string f="exp(x)-2*x-1";
    ex x0=0;
    ex x1=1;
    ex x2=2;
    ex tol=1e-10;
    int iterMax=2500;
    vector<ex> valores=muller(f,x0,x1,x2,tol,iterMax);
    cout<<"El valor aproximado de x es "<<valores[0].evalf()<<" con un error de "<<valores[1]<<" encontrado en "<<valores[2]<<" iteraciones"<<endl;

    return 0;
}