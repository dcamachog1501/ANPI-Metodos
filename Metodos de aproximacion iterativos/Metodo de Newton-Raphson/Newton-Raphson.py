def newton_raphson(f,x0,tol,iterMax):
    x=symbols('x')
    f=sympify(f)
    df=diff(f,x)
    if(df!=0):
        f=lambdify(x,f)
        df=lambdify(x,df)
        k=0
        xk=x0
        error=None
        while(k<iterMax):
            error=abs(f(xk))
            if(error<=tol):
                return(xk,error,k)
            if(df(xk)==0):
                print("Indefinicion en alguna de las iteraciones, se retornan los valores hasta dicha iteracion")
                return(xk,error,k)
            xk=xk-(f(xk)/df(xk))
            k+=1
        return(xk,error,k)
    print("Esta funcion no puede ser aproximada por el metodo de Newton-Raphson,la derivada de f(x) es 0")
    xk=None
    return(xk,error,k)

f="cos(2*x)**2-x**2"
valores=newton_raphson(f,100,10**-100,2500)
if(valores[0]!=None):
    print("La aproximacion encontrada es x= "+str(valores[0])+" con un error de "+str(valores[1])+" en "+str(valores[2])+" iteraciones")
