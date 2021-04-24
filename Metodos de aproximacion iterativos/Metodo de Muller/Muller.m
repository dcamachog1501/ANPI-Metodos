function Muller
  pkg load symbolic
  tol=1e-10;
  iterMax=2500;
  x0=0;
  x1=1;
  x2=2;
  f="exp(x)-2*x-1";
  [xk,error,k]=muller(f,x0,x1,x2,tol,iterMax);
  if(xk!=inf)
    display(["El valor aproximado por el metodo es de ",num2str(xk)," con un error de ",num2str(error)," en ",num2str(k)," iteraciones"]);
  endif
end

function[xk,error,k]=muller(f,x0,x1,x2,tol,iterMax)
  f=sym(f);
  f=matlabFunction(f);
  xk=inf;
  error=inf;
  k=0;
  for i=0:iterMax
    
    if(x0-x1==0||x0-x2==0||x1-x2==0)
      disp("Indefinicion en alguna de las iteraciones, se retornan los valores hasta dicha iteracion");
      return
    endif
    a = ((x1 - x2) * (f(x0) - f(x2)) - (x0 - x2) * (f(x1) - f(x2))) / ((x0 - x1) * (x0 - x2) * (x1 - x2));
    b = (((x0 - x2) ** 2) * (f(x1) - f(x2)) - ((x1 - x2) ** 2) * (f(x0) - f(x2))) / ((x0 - x1) * (x0 - x2) * (x1 - x2));
    c = f(x2);
    xk=x2-((2*c)/(b+(b/abs(b))*sqrt(b**2-4*a*c)));
    error=f(xk);
    
    if(abs(error)<=tol)
      error=abs(error);
      return
    endif
    
    if(abs(x0-xk)<abs(x2-xk))
      x2=x1;
    else
      x0=x1;
    endif
    x1=xk;
    k=k+1;
  endfor
  error=abs(error);
  return
end