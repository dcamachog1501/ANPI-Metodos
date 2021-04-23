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
ex funct_aux(string f,ex x_0,ex x_1)
{
    ex f_aux=x_1-(((x_1-x_0)/(funct(f,x_1)-funct(f,x_0)))* funct(f,x_1));
    return evalf(f_aux);
}
vector<ex>secante(string f,ex x0,ex x1,ex tol,int iterMax)
{
    vector<ex> output;
    ex xk=x1;
    ex xk_1=x0;
    ex error=0;
    int k=0;
    while(k<iterMax)
    {
        if(funct(f,xk_1)- funct(f,xk)==0)
        {
            cout<<"Indefinicion en alguna de las iteraciones,se retornan los valores hasta dicha iteracion"<<endl;
            break;
        }
        error= funct(f,xk);
        if(abs(error)<=tol)
        {
            break;
        }
        ex xkaux= funct_aux(f,xk_1,xk);
        xk_1=xk;
        xk=xkaux;
        k++;
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
    ex x1=200;
    ex tol=1e-100;
    int iterMax=2500;
    vector<ex> valores=secante(f,x0,x1,tol,iterMax);
    if(valores[0]!='x')
    {
        cout<<"El valor aproximado de x es "<<valores[0].evalf()<<" con un error de "<<valores[1]<<" encontrado en "<<valores[2]<<" iteraciones"<<endl;
    }
    return 0;
}