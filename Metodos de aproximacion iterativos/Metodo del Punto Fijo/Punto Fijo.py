import numpy as np
from sympy import *
from math import *

def punto_fijo(f,g,x0,tol,iterMax):
    x=symbols('x')

    f=lambdify(x,f)
    g=lambdify(x,g)
    xk=x0
    error=0
    k=0
    while(k<iterMax):
        error=f(xk)
        if(abs(error)<=tol):
            return(xk,abs(error),k)
        k+=1
        xk = g(xk)
    return(xk,abs(error),k)

f="cos(2*x)**2-x**2"
g="cos(2*x)"

valores=punto_fijo(f,g,0,1e-100,25000)

if(valores[0]!=None):
    print("La aproximacion encontrada es x= "+str(valores[0])+" con un error de "+str(valores[1])+" en "+str(valores[2])+" iteraciones")
