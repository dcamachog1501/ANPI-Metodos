import numpy as np
from math import *

def sust_atras(a,b):
    n = len(b)
    xk=np.zeros((1,n))
    for i in range(n-1,-1,-1):
        suma=0
        for j in range(i+1,n):
            suma+=a[i,j]*xk[0,j]
        xi=(1/a[i,i])*(b[i]-suma)
        xk[0,i]=xi
    return xk
def sust_adelante(a,b):
    n=len(b)
    xk = np.zeros((1, n))
    for i in range(n):
        suma=0
        for j in range(i):
            suma+=a[i,j]*xk[0,j]
        xi=(1/a[i,i])*(b[i]-suma)
        xk[0,i]=xi
    return xk

def is_symetric(a):
    dims=a.shape
    a_t=a.transpose()
    for i in range(dims[0]):
        for j in range(dims[1]):
            if(a[i,j]!=a_t[i,j]):
                return False
    return True
def cholesky_l(a):
    dims = a.shape
    n=dims[0]
    if (dims[0]==dims[1]):
        for i in range(n):
            temp = a[0:i+1, 0:i+1]
            if (np.linalg.det(temp)<0):
                print("Uno de los deteminates de las submatrices principales de a no es positivo")
                return None
        if(not is_symetric(a)):
            print("La matriz provista no es simetrica")
            return None

        l=np.zeros(dims)
        for i in range(n):
            for j in range(i+1):
                sum=0;
                #Calculo del valor de la diagonal
                if(j==i):
                    for k in range(j):
                        sum+=pow(l[j,k],2)
                    l[j,j]=sqrt(a[j,j]-sum);
                #Calculo de los valores debajo de la diagonal
                else:
                    for k in range(j):
                        sum+=l[i,k]*l[j,k];
                    l[i,j]=(a[i,j]-sum)/l[j,j]
        return l
    print("La matriz no es cuadrada")
    return None

def fact_cholesky(a,b):
    l=cholesky_l(a)
    if(l is not None):
        yk=sust_adelante(l,b)
        return sust_atras(l.transpose(),yk.transpose())
    at=np.matmul(a.transpose(),a)
    bt=np.matmul(a.transpose(),b)
    print("Aplicando substitucion por A techo...\n")
    return fact_cholesky(at,bt)

a=np.matrix(([1,-1,2],
             [-2,0,4],
             [0,-2,7]))
b=np.matrix(([0],
             [2],
             [5]))

xk=fact_cholesky(a,b)
if(xk is not None):
    print("Las soluciones del sistema de ecuaciones son: \nxk="+str(xk))
