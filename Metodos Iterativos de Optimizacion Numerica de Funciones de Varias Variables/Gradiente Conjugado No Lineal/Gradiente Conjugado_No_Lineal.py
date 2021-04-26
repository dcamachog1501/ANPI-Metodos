import numpy as np
from sympy import *
from math import *
import random

def gradient(f,vars):
    output=[]
    for var in vars:
        output.append(sympify(diff(f,var)))
    return output
def gradientEval(gradient,vars,values):
    output=[]
    for ind in range(len(vars)):
        output.append(gradient[ind].subs(zip(vars,values)))
    return output
def norm(gradient):
    suma=0
    for component in gradient:
        suma+=component**2
    return sqrt(suma)
def gcnl(f,variables,iniciales,tol,iterMax):
    var_temp=variables.split(" ")
    val_temp=iniciales.split(" ")
    vars=[]
    xk=[]
    for ind in range(len(var_temp)):
        vars.append(var_temp[ind])
        xk.append(float(val_temp[ind]))

    f=sympify(f)
    grdf=gradient(f,vars)
    gk=gradientEval(grdf,vars,xk)
    dk=[]
    for valor in gk:dk.append(-valor)
    k,error=0,0
    while(k<iterMax):
        delta=random.random()
        ak=0.5
        aux=[]
        for ind in range(len(xk)):aux.append(xk[ind]+ak*dk[ind])
        izq=f.subs(zip(vars,aux))-f.subs(zip(vars,xk))
        der=0
        for ind in range(len(gk)):der+=delta*ak*gk[ind]*dk[ind]
        while(izq>der):
            ak=ak/2
            aux = []
            for ind in range(len(xk)): aux.append(xk[ind] + ak * dk[ind])
            izq = f.subs(zip(vars, aux)) - f.subs(zip(vars, xk))
            der = 0
            for ind in range(len(gk)): der += delta * ak * gk[ind] * dk[ind]
        xkp1=[]
        for ind in range(len(xk)):xkp1.append(xk[ind]+ak*dk[ind])
        xk=xkp1
        error=norm(gradientEval(grdf,vars,xk))
        if(abs(error)<=tol):
            return (xk,abs(error),k)
        gkp1=gradientEval(grdf,vars,xk)
        beta=(((norm(gkp1))**2)/((norm(gk))**2))
        gk=gkp1
        dkp1=[]
        for ind in range(len(dk)):dkp1.append(-gk[ind]+beta*dk[ind])
        dk=dkp1
        k+=1
    return (xk, abs(error), k)


f="(x-2)^4+(x-2*y)^2"
variables="x y"
iniciales="3 0"#Componentes del vector inicial
tol=1e-5
iterMax=2500
valores=gcnl(f,variables,iniciales,tol,iterMax)
print("El punto aproximado es "+str(valores[0])+" con un error de "+str(valores[1])+" en "+str(valores[2])+" iteraciones")
