import numpy as np
from math import *


def triang_sup(a,b):
    n=len(b)
    a_aumt=np.concatenate((a,b),axis=1)
    for k in range(0,n-1):
        for i in range(k+1,n):
            mik=a_aumt[i,k]/a_aumt[k,k]
            for j in range(k,n+1):
                a_aumt[i,j]=a_aumt[i,j]-mik*a_aumt[k,j];
                print(mik)
                print(a_aumt)
    at=a_aumt[0:n+1,0:n]
    bt=a_aumt[0:n+1,n]

    return (at,bt)
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

def gaussiana(a,b):
    if(np.linalg.det(a)!=0):
        triang=triang_sup(a,b)
        return sust_atras(triang[0],triang[1])
    return ([None])
    

soluciones=gaussiana(np.matrix(([2,-6,12,16],[1,-2,6,6],[-1,3,-3,-7],[0,4,3,-6])),np.matrix('70;26;-30;-26'))
if(soluciones[0,0]!=None):
    i=1
    sol=""
    for solucion in soluciones:
        sol+="x"+str(i)+"="+str(solucion)+" "
    print("Soluciones encontradas para el sistema de ecucaciones:\n"+sol)
