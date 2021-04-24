import numpy as np
from sympy import *
from math import *

def muller(f,x0,x1,x2,tol,iterMax):
    x=symbols('x')
    f=lambdify(x,f)
    xk=0
    error=0
    k=0
    while(k<iterMax):
        if(x0-x1==0 or x0-x2==0 or x1-x2==0):
            print("Indefinicion en alguna de las iteraciones, se retornan los valores hasta dicha iteracion")
            break
        a = ((x1 - x2) * (f(x0) - f(x2)) - (x0 - x2) * (f(x1) - f(x2))) / ((x0 - x1) * (x0 - x2) * (x1 - x2))
        b = (((x0 - x2) ** 2) * (f(x1) - f(x2)) - ((x1 - x2) ** 2) * (f(x0) - f(x2))) / ((x0 - x1) * (x0 - x2) * (x1 - x2))
        c = f(x2)
        xk=x2-((2*c)/(b+(b/abs(b))*sqrt(b**2-4*a*c)))
        error=f(xk)
        if(abs(error)<=tol):
            break
        if(abs(xk-x0)<abs(xk-x1)):
            x2=x1
        else:
            x0=x1

        x1=xk
        k+=1
    return (xk,abs(error),k)

f="exp(x)-2*x-1"
valores=muller(f,0,1,2,1e-10,2500)
if(valores[0]!=None):
    print("La aproximacion encontrada es x= "+str(valores[0])+" con un error de "+str(valores[1])+" en "+str(valores[2])+" iteraciones")
