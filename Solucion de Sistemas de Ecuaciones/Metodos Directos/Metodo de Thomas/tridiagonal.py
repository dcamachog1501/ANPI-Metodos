import numpy as np

def tridiagonal(a,b,c,n):
    out=np.zeros((n,n))
    for i in range(n):
        if(i==0):
            out[i,i+1]=c
        elif(i==n-1):
            out[i,i-1]=a
        else:
            out[i,i+1]=c
            out[i,i-1]=a
        out[i,i]=b
    return out
