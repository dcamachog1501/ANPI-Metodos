import numpy as np

def tridiagonal(a,b,c,n):
    out=np.zeros((n,n))
    for i in range(n):
        if(i==0):
            out[i,i+1]=b
        elif(i==n-1):
            out[i,i-1]=c
        else:
            out[i,i+1]=b
            out[i,i-1]=c
        out[i,i]=a
    return out
