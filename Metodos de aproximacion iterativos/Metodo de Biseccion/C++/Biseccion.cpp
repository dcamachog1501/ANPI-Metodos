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
vector<ex>biseccion(ex a,ex b,string f,ex tol,int iterMax)
{
    vector<ex> output;
    if(funct(f,a)*funct(f,b)<0)
    {
        int k=0;
        ex xk;
        ex error;
        while(k<iterMax)
        {
            xk=evalf((a+b)/2);
            error=funct(f,xk);
            if(abs(error)<=tol)
            {
                break;
            }
            if(funct(f,a)*error<0)
            {
                b=xk;
            }
            else if(funct(f,b)*error<0)
            {
                a=xk;
            }
            else
            {
                cout<<"Ninguno de los nuevos intervalos cumple con el teorema de Bolzano"<<endl;
                output.push_back('x');
                return output;
            }
            k++;
        }
        output.push_back(xk);
        output.push_back(abs(error));
        output.push_back(k);
        return output;
    }
    cout<<"La funcion no cumple con el teorema de Bolzano en el intervalo especificado"<<endl;
    output.push_back('x');
    return output;
}

int main()
{
    Digits=10;
    string f="pow(cos(2*x),2)-pow(x,2)";
    ex a=0;
    ex b=Pi;
    ex tol=1e-100;
    int iterMax=2500;
    vector<ex> valores=biseccion(a,b,f,tol,iterMax);
    if(valores[0]!='x')
    {
        cout<<"El valor aproximado de x es "<<valores[0].evalf()<<" con un error de "<<valores[1]<<" encontrado en "<<valores[2]<<" iteraciones"<<endl;
    }
    return 0;
}
