import numpy as np
from math import *

def fact_lu(a,b):
    lu_mat=lu(a)
    if lu_mat!=None:
        yk=sust_adelante(lu_mat[0],b)
        return sust_atras(lu_mat[1],yk.transpose())
    return None
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
def lu(a):
    dims=a.shape
    n=len(a)
    if(dims[0]==dims[1]):
        for i in range(n):
            temp=a[0:i+1,0:i+1]
            if(np.linalg.det(temp)==0):
                print("Una de las submatrices de a no es invertible")
                return None
        l=np.identity(n)
        for k in range(n):
            for i in range(k+1,n):
                l[i,k]=a[i,k]/a[k,k]
                a[i,k]=0
                for j in range(k+1,n):
                    a[i,j]=a[i,j]-l[i,k]*a[k,j];
        u=a
        return(l,u)
    else:
        print("La matriz provista no es cuadrada")
        return None

a=np.matrix(([4 ,-2, 1],
             [20,-7,12],
             [-8,13,17]))
b=np.matrix(([11],
             [70],
             [17]))
xk=fact_lu(a,b)
if(xk is not None):
    print("Las soluciones dle sitema de ecuaciones son: \nxk="+str(xk))
