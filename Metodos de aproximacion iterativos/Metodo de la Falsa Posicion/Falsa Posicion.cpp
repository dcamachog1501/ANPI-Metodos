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
    ex f_aux=x_1-(((x_1-x_0)/(funct(f,x_1)-funct(f,x_0))))* funct(f,x_1);
    return evalf(f_aux);
}
vector<ex>falsa_posicion(string f,ex a,ex b,ex tol,int iterMax)
{
    vector<ex> output;
    if(funct(f,a)* funct(f,b)<0)
    {
        ex xk = 0;
        ex error = 0;
        int k = 0;
        while (k < iterMax)
        {
            if (funct(f, b) - funct(f, a) == 0)
            {
                cout << "Indefinicion en alguna de las iteraciones,se retornan los valores hasta dicha iteracion"<< endl;
                break;
            }
            xk = funct_aux(f, a, b);
            error = funct(f, xk);
            if (abs(error)<= tol)
            {
                break;
            }
            if (funct(f, a) * error < 0)
            {
                b = xk;
            } else if (funct(f, b) * error < 0)
            {
                a = xk;
            }
            else
            {
                cout << "Ninguno de los nuevos intervalos cumple con el teorema de Bolzano" << endl;
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
    output.push_back('x');
    return output;
}

int main()
{
    Digits=10;
    string f="pow(cos(2*x),2)-pow(x,2)";
    ex a=0;
    ex b=Pi;
    ex tol=1e-10;
    int iterMax=2500;
    vector<ex> valores=falsa_posicion(f,a,b,tol,iterMax);
    if(valores[0]!='x')
    {
        cout<<"El valor aproximado de x es "<<valores[0].evalf()<<" con un error de "<<valores[1]<<" encontrado en "<<valores[2]<<" iteraciones"<<endl;
    }
    return 0;
}