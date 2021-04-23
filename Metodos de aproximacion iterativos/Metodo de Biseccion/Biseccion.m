function Biseccion
  f='(cos(2*x)^2)-x^2';
  a=0;
  b=pi;
  x0=100;
  tol=1*10^-100;
  iterMax=2500;
  [xk,error,k]=biseccion_aux(f,a,b,tol,iterMax);
  if(xk!=inf)
    display(["El valor aproximado por el metodo es de ",num2str(xk)," con un error de ",num2str(error)," en ",num2str(k)," iteraciones"])
  endif
end
  
function[xk,error,k]=biseccion_aux(f,a,b,tol,iterMax)
  pkg load symbolic
  
  f=sym(f);
  f=matlabFunction(f);
  xk=inf;
  error=inf;
  k=0;
  if(f(a)*f(b)<0)
    for i=1:iterMax
      xk=(a+b)/2;
      error=f(xk);
      if(abs(error)<=tol)
      error=abs(error);
        return
      endif
      
      if(error*f(a)<0)
        b=xk;
        
      elseif(error*f(b)<0)
        a=xk;
        
      else
        disp("Ninguno de los nuevos intervalos cumple con el teorema de Bolzano");
        xk=inf;
        return
      endif
      k=k+1;
    endfor
    error=abs(error);
    return
  endif
  disp("La funcion no puede ser aproximada por el metodo de la biseccion, no cumple con el teorema de Bolzano")
  return
 endfunction