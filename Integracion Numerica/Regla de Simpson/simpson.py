from sympy import *
import numpy as np
import scipy.optimize as opt

def simpson(f,a,b):
    x=symbols('x')
    h=(b-a)/2

    d_4=diff(f,n=4)

    #Funcion para encontrar el maximo en [a,b] de d_4(f_4=-|d_4| minimo de f_4= maximo de |d_4|)
    f_4=-1*abs(d_4)
    f_4=lambdify(x,f_4)

    f=lambdify(x,f)
    d_4=lambdify(x,d_4)

    I=(f(a)+4*f((a+b)/2)+f(b))*h/3

    beta=abs(d_4(opt.minimize_scalar(f_4,bounds=[a,b],method='bounded').x))

    error=beta*((h**5)/90)

    return I,error

f='ln(x)'
a=2
b=5
I,error=simpson(f,a,b)

print("El valor aproximado de la integral es  "+str(I)+" con un error maximo de "+str(error))
