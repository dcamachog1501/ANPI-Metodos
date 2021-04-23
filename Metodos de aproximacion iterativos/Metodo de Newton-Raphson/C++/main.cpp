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
ex dfunct(string f,ex value)
{
    symbol x;
    symtab table;
    table["x"]=x;
    parser reader(table);
    ex expresion=reader(f);
    expresion=expresion.diff(x,1);
    return evalf(expresion.subs(x==value));
}
vector<ex>newton_raphson(string f,ex x0,ex tol,int iterMax)
{
    vector<ex> output;
    ex xk=x0;
    ex error=0;
    int k=0;
    while(k<iterMax)
    {
        error= funct(f,xk);
        if(abs(error)<=tol)
        {
            break;
        }
        k++;
        ex diffxk=dfunct(f,xk);
        if(diffxk==0)
        {
            cout<<"Indefinicion en alguna de las iteraciones,se retornan los valores hasta dicha iteracion"<<endl;
            break;
        }
        xk=xk-(funct(f,xk)/dfunct(f,xk));
    }
    output.push_back(xk);
    output.push_back(error);
    output.push_back(k);
    return output;
}

int main()
{
    Digits=10;
    string f="pow(cos(2*x),2)-pow(x,2)";
    ex x0=100;
    ex tol=1e-10;
    int iterMax=2500;
    vector<ex> valores=newton_raphson(f,x0,tol,iterMax);
    if(valores[0]!='x')
    {
        cout<<"El valor aproximado de x es "<<valores[0].evalf()<<" con un error de "<<valores[1]<<" encontrado en "<<valores[2]<<" iteraciones"<<endl;
    }
    return 0;
}