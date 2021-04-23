import numpy as np
from sympy import *
from math import *

from sympy.calculus.util import continuous_domain

#Metodo para determinar si una funcion es contiua en el intervalo de [a,b]
'''Entradas: f=Funcion a evaluar
             a=Limite inferior del intervalo en el que se debe evaluar la continuidad de la funcion
             b=limite superior del intervalo en el que se debe evaluar la continuidad de la funcion
   Salidas: Booleano que inidica si la funcion es continua sobre el intervalo [a,b]'''

def isContinuous(f,a,b):
    f=sympify(f)
    x=symbols('x')
    var=continuous_domain(f,x,Interval(a,b))
    return Interval(a,b)==var


def biseccion(a,b,f,tol,iterMax):
    if(isContinuous(f,a,b)):
        x=Symbol('x')
        f=lambdify(x,f)
        if(f(a)*f(b)<0):
            k=0
            xk=0
            error=0
            while(k<=iterMax):
                xk=(a+b)/2
                error=f(xk)
                if(abs(error)<=tol):
                    return (xk,abs(error),k)
                if(f(a)*error<0):
                    b=xk
                elif(f(b)*error<0):
                    a=xk
                else:
                    print("Ninguno de los nuevos intervalos cumple con el teorema de Bolzano")
                    return (None,None,k)
                k+=1
            return (xk,abs(error),k)
        else:
            print("No es posible utilizar el metodo de biseccion para esta funcion, puesto que no cumple con el teorema de Bolzano")
            return (None,None,None)

    print("La funcion no puede ser evaluada puesto que no es continua")
    return (None,None,None)

f="cos(2*x)**2-x**2"
valores=biseccion(0,pi,f,10**-100,2500)
if(valores[0]!=None):
    print("La aproximacion encontrada es x= "+str(valores[0])+" con un error de "+str(valores[1])+" en "+str(valores[2])+" iteraciones")
