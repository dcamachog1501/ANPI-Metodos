from sympy import *
import numpy as np
import scipy.optimize as opt



def simpson_compuesto(f,a,b,n):
    if(n//2!=0):
        vx=np.linspace(a,b,n)
        I=0

        for i in range(n-1):
            I+=simpson(f,vx[i],vx[i+1])

        h=(b-a)/(n-1)

        x=symbols('x')
        d_4 = diff(f, n=4)
        f_4 = -1 * abs(d_4)

        f_4 = lambdify(x, f_4)
        d_4 = lambdify(x, d_4)

        beta = abs(d_4(opt.minimize_scalar(f_4, bounds=[a, b], method='bounded').x))

        error=beta*((b-a)*h**4)/180
        return (I,error)

    print("No se puede trabajar con un numero n par de puntos")
    return None

def simpson(f, a, b):
    x = symbols('x')
    h = (b - a) / 2

    f = lambdify(x, f)

    return (f(a) + 4 * f((a + b) / 2) + f(b)) * h / 3

f='ln(x)'

a=2
b=5
n=239

values=simpson_compuesto(f,a,b,n)
if(values!=None):
    print("El valor aproximado es de "+str(values[0])+" con un error maximo de "+str(values[1]))
