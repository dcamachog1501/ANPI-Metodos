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
vector<ex>punto_fijo(string f,string g,ex x0,ex tol,int iterMax)
{
    vector<ex> output;
    ex xk = x0;
    ex error = 0;
    int k = 0;
    while (k < iterMax)
    {
        error=funct(f,xk);
        if (abs(error)<= tol)
        {
            break;
        }
        xk=funct(g,xk);
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
    string g="log(2*x+1)";
    ex x0=1.5;
    ex tol=1e-10;
    int iterMax=2500;
    vector<ex> valores=punto_fijo(f,g,x0,tol,iterMax);
    cout<<"El valor aproximado de x es "<<valores[0].evalf()<<" con un error de "<<valores[1]<<" encontrado en "<<valores[2]<<" iteraciones"<<endl;

    return 0;
}