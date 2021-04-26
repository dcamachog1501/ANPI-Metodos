import numpy as np
import scipy.optimize as opt
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
def dc(f,variables,iniciales,tol,iterMax):
    var_temp=variables.split(" ")
    val_temp=iniciales.split(" ")
    vars=[]
    xk=[]
    for ind in range(len(var_temp)):
        vars.append(var_temp[ind])
        xk.append(float(val_temp[ind]))

    f=sympify(f)
    grdf=gradient(f,vars)
    error=0
    for k in range (iterMax):
        for j in range(len(vars)):
            var_temp=[]
            val_temp=[]
            for i in range(len(vars)):
                if(i!=j):
                    var_temp.append(vars[i])
                    val_temp.append(xk[i])
            f_temp=f.subs(zip(var_temp,val_temp))
            xk[j]=opt.fmin(lambdify(vars[j],f_temp),0,disp=False)[0]
        error=norm(gradientEval(grdf,vars,xk))
        print(abs(error))
        if(error<=tol):
            return (xk,error,k)
    return (xk,error,k)


f="(x-2)^4+(x-2*y)^2"
variables="x y"
iniciales="-0.5 -0.5"#Componentes del vector inicial
tol=1e-2
iterMax=2500
valores=dc(f,variables,iniciales,tol,iterMax)
print("El punto aproximado es "+str(valores[0])+" con un error de "+str(valores[1])+" en "+str(valores[2])+" iteraciones")