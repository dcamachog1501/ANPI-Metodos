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

def falsa_posicion(f,a,b,tol,iterMax):
    if(isContinuous(f,a,b)):
        x=symbols('x')
        f=lambdify(x,f)
        if(f(a)*f(b)<0):
            xk,error,k=0,0,0
            while(k<iterMax):
                if(f(b)-f(a)==0):
                    print("Indefinicion en alguna de las iteraciones, se retornan los valores hasta dicha iteracion")
                    break
                xk = b-((b - a)/(f(b) - f(a)))* f(b)
                error=f(xk)
                if(abs(error)<=tol):
                    break
                if(error*f(a)<0):
                    b=xk
                elif(error*f(b)<0):
                    a=xk
                else:
                    print("Ninguno de los nuevos intervalos cumple con el teorema de Bolzano")
                    return (None,None,None)
                k+=1
            return(xk,error,k)
        else:
            print("La funcion no cumple con el teorema de Bolzano")
            return (None,None,None)
    else:
        print("La funcion provista no es continua en el intervalo espcificado")
        return (None,None,None)
f="sin(x)**2+x**2-1"
valores=falsa_posicion(f,0,2,10**-10,2500)
if(valores[0]!=None):
    print("La aproximacion encontrada es x= "+str(valores[0])+" con un error de "+str(valores[1])+" en "+str(valores[2])+" iteraciones")
