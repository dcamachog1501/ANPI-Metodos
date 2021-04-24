import numpy as np
from sympy import *
from math import *

def secante(f,x0,x1,tol,iterMax):
    x = symbols('x')
    f=lambdify(x,f)
    xk=x1
    k=0
    error=0
    while(k<=iterMax):
        error=f(xk)
        if(abs(error)<=tol):
            break
        if(f(x1)-f(x0)==0):
            print("Indefinicion en alguna de las iteraciones, se retornan los valores hasta dicha iteracion")
            break
        xk=x1-((x1-x0)/(f(x1)-f(x0)))*f(x1)
        x0=x1
        x1=xk
        k+=1
    return (xk, error, k)

f="cos(2*x)**2-x**2"
valores=secante(f,100,200,10**-100,2500)
if(valores[0]!=None):
    print("La aproximacion encontrada es x= "+str(valores[0])+" con un error de "+str(valores[1])+" en "+str(valores[2])+" iteraciones")
