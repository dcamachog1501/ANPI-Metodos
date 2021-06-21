from numpy import *
from sympy import *
import scipy.optimize as opt

def trapecio(f,a,b):
    x=symbols('x')
    f=lambdify(x,f)
    h=b-a
    I=h*(f(a)+f(b))/2
    return I

def trapecio_compuesto(f,a,b,n):
    x = symbols('x')
    xv=linspace(a,b,n)#Funcion para crear un vector de n valores entre a y b equidistantes(no es necesario el h)--->Equivalente en Octave: xv=a:h:b
    I=0
    for i in range(n-1):
        I+=trapecio(f,xv[i],xv[i+1])


    f = sympify(f)
    df_2 = diff(diff(f, x), x)

    f_2 = lambdify(x, -abs(df_2))
    df_2 = lambdify(x, df_2)
    alpha = abs(df_2(opt.minimize_scalar(f_2, bounds=[a,b], method='bounded').x))
    h=(b-a)/(n-1)

    error=alpha*((b-a)*h**2)/12


    return (I,error)

f='ln(x)'
values=trapecio_compuesto(f,2,5,239)
print("El valor aproximado es de "+str(values[0])+" con un error maximo de "+str(values[1]))
