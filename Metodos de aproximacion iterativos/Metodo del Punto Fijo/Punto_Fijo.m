function Punto_Fijo
  pkg load symbolic
  tol=1e-10;
  iterMax=2500;
  x0=1.5;
  f="exp(x)-2*x-1"
  g="ln(2*x+1)"
  [xk,error,k]=pfijo(f,g,x0,tol,iterMax);
  if(xk!=inf)
    display(["El valor aproximado por el metodo es de ",num2str(xk)," con un error de ",num2str(error)," en ",num2str(k)," iteraciones"])
  endif
end

function[xk,error,k]=pfijo(f,g,x0,tol,iterMax)
  f=sym(f);
  g=sym(g);
  f=matlabFunction(f);
  g=matlabFunction(g);
  xk=x0;
  error=inf;
  k=0;
  
  for i=0:iterMax
    error=f(xk);
    if(abs(error)<=tol)
      return
    endif
    xk=g(xk);
    k=k+1;
  endfor
  return

end