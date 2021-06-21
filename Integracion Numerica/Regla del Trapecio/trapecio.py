import numpy as np
import scipy as scp
from sympy import *
import scipy.optimize as opt

def trapecio(f,a,b):
    x= symbols('x')
    f=sympify(f)
    df_2 =diff(diff(f,x),x)

    f=lambdify(x,f)
    f_2 = lambdify(x, -abs(df_2))
    df_2=lambdify(x,df_2)


    h=b-a
    I=h*(f(a)+f(b))/2


    alpha=abs(df_2(opt.minimize_scalar(f_2,bounds=[a,b],method='bounded').x))
    error=alpha*(h**3)/12

    return (I,error)

f='ln(x)'
a=2
b=5
values=trapecio(f,a,b)

print("El valor aproximado de la integral es  "+str(values[0])+" con un error maximo de "+str(values[1]))