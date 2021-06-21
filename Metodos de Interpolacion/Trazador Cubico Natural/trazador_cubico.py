import numpy as np
from sympy import *

def trazadorCubico(points):
    # Generacion del vector de intervalos h
    h=[]
    for ind in range(len(points)-1):
        h.append(points[ind+1][0]-points[ind][0])
    xs=[]
    ys=[]
    for point in points:
        xs.append(point[0])
        ys.append(point[1])

    # Generacion de la matriz tridiagonal A
    A=tridiagonal(h,len(points)-2)

    # Generacion del vector u
    u=np.zeros((1,len(points)-2))
    for i in range(len(points)-2):
        u[0][i]=(6*((ys[i+2]-ys[i+1])/(h[i+1])-(ys[i+1]-ys[i])/(h[i])))
    u=u.transpose()

    # Resolviendo A*M=u y agregando M0 y Mn
    M=thomas(A,u)
    M=np.append(np.append([0],M),[0])

    trazadores = []
    x = symbols('x')

    # Se calculan los coeficientes a,b,c,d y se calculan los traadores Si
    for i in range(len(points)-1):
        a=(M[i+1]-M[i])/6*h[i]
        b=M[i]/2
        c=((ys[i+1]-ys[i])/h[i])-((h[i]*(M[i+1]+2*M[i]))/6)
        d=ys[i]

        trazadores.append(a*(x-xs[i])**3+b*(x-xs[i])**2+c*(x-xs[i])+d)

    return trazadores

def tridiagonal(h,n):
    out=np.zeros((n,n))
    for i in range(n):
        if(i==0):
            out[i,i+1]=h[i+1]
        elif(i==n-1):
            out[i,i-1]=h[i]
        else:
            out[i,i+1]=h[i+1]
            out[i,i-1]=h[i]
        out[i,i]=2*(h[i]+h[i+1])
    return out

def is_tridiagonal(a):
    dims=a.shape

    for i in range(dims[0]):
        for j in range(dims[1]):
            if(j>i+1 and a[i,j]!=0):
                return False
            elif(j<i-1 and a[i,j]!=0):
                return False
    return True
def thomas(a,b):
    if(is_tridiagonal(a)):
        n=len(b)
        a_s=[]
        b_s=[]
        c_s=[]
        d_s=b.transpose();
        p_s=[]
        q_s=[]
        xk =np.zeros((1,n))
        for i in range(n):
            if(i==0):
                a_s.append(0)
                b_s.append(a[i,i])
                c_s.append(a[i,i+1])
            elif(i==n-1):
                a_s.append(a[i,i-1])
                b_s.append(a[i,i])
                c_s.append(0)
            else:
                a_s.append(a[i, i - 1])
                b_s.append(a[i, i])
                c_s.append(a[i, i + 1])
        n=len(b)
        qi=0
        pi=0
        for i in range(n):
            if(i==0):
                p_s.append(c_s[i]/b_s[i])
                q_s.append(d_s[0,i]/b_s[i])

            else:
                if(i!=n-1):
                    p_s.append(c_s[i]/(b_s[i]-p_s[i-1]*a_s[i]))
                    q_s.append((d_s[0,i]-q_s[i-1]*a_s[i])/(b_s[i]-p_s[i-1]*a_s[i]))
                else:
                    q_s.append((d_s[0, i] - q_s[i - 1] * a_s[i]) / (b_s[i] - p_s[i - 1] * a_s[i]))
        for j in range(n-1,-1,-1):
            if(j==n-1):
                xk[0,j]=q_s[j]
            else:
                xk[0,j]=q_s[j]-p_s[j]*xk[0,j+1]
        return xk
    else:
        print("La matriz no es tridiagonal")
        return None

points=[(1,2.718282),(1.05,3.286299),(1.07,3.527609),(1.1,3.905416)]
print(trazadorCubico(points))