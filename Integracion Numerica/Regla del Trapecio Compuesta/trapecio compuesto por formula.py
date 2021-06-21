from numpy import *
from sympy import *
import scipy.optimize as opt

def trapecio_compuesto(f,a,b,n):
    x = symbols('x')
    xv=linspace(a,b,n)
    h=(b-a)/(n-1)

    f = sympify(f)
    df_2 = diff(diff(f, x), x)

    f_2 = lambdify(x, -abs(df_2))
    df_2 = lambdify(x, df_2)
    f=lambdify(x,f)


    I=0
    for i in range(n-1):
        I+=f(xv[i])+f(xv[i+1])
    I=I*h/2

    alpha = abs(df_2(opt.minimize_scalar(f_2, bounds=[a,b], method='bounded').x))
    h=(b-a)/(n-1)

    error=alpha*((b-a)*h**2)/12


    return (I,error)

f='ln(x)'
values=trapecio_compuesto(f,2,5,239)
print("El valor aproximado es de "+str(values[0])+" con un error maximo de "+str(values[1]))
