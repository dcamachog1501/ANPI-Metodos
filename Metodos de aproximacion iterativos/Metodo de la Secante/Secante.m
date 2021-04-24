function Secante
  pkg load symbolic;
  f='(cos(2*x)^2)-x^2';
  a=0;
  b=pi;
  x0=100;
  x1=200;
  tol=1*10^-100;
  iterMax=2500;
  [xk,error,k]=secante_aux(f,x0,x1,tol,iterMax);
  if(xk!=inf)
    display(["El valor aproximado por el metodo es de ",num2str(xk)," con un error de ",num2str(error)," en ",num2str(k)," iteraciones"])
  endif
end

function[xk,error,k]=secante_aux(f,x0,x1,tol,iterMax)
  f=matlabFunction(sym(f));
  xk=x1;
  error=0;
  k=0;
  for i=1:iterMax
    error=f(xk);
    if(abs(error)<=tol)
      return
    endif
    if(f(x1)-f(x0)==0)
      disp("Indefinicion en alguna de las iteraciones, se retornan los valores hasta dicha iteracion");
      return
    endif
    xk=x1-((x1-x0)/(f(x1)-f(x0)))*f(x1);
    x0=x1;
    x1=xk;
    k=k+1;
  endfor
  return
end