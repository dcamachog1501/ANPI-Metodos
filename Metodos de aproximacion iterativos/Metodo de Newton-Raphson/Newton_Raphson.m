function Newton_Raphson
  pkg load symbolic
  f='(cos(2*x)^2)-x^2';
  x0=100;
  tol=1*10^-100;
  iterMax=2500;
  [xk,error,k]=newton_raphson(f,x0,tol,iterMax);
  if(xk!=inf)
    display(["El valor aproximado por el metodo es de ",num2str(xk)," con un error de ",num2str(error)," en ",num2str(k)," iteraciones"])
  endif
end
  
function[xk,error,k]=newton_raphson(f,x0,tol,iterMax)
  f=sym(f);
  x=sym('x');
  df=diff(f,x);
  if(df!=0)
    f=matlabFunction(f);
    df=matlabFunction(df);
    
    xk=x0;
    error=inf;
    k=0;
    for i=1:iterMax
      error=abs(f(xk));
      if(error<=tol)
        return
      endif
      
      xk=xk-(f(xk)/df(xk));
      k=k+1;
    endfor
  endif
  disp('No es posible aproximar esta funcion con el metodo de Newton-Raphson,su derivada es 0')
  return
  
end