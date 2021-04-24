function Falsa_Posicion
  pkg load symbolic
  a=0;
  b=pi;
  tol=1e-10;
  iterMax=2500;
  f="cos(2*x)^2-x^2";
  [xk,error,k]=falsa(f,a,b,tol,iterMax);
  if(xk!=inf)
    display(["El valor aproximado por el metodo es de ",num2str(xk)," con un error de ",num2str(error)," en ",num2str(k)," iteraciones"])
  endif
end

function[xk,error,k]=falsa(f,a,b,tol,iterMax)
  f=sym(f);
  f=matlabFunction(f);
  xk=inf;
  error=inf;
  k=0;
  if(f(a)*f(b)<0)
    for i=0:iterMax
      if(f(b)-f(a)==0)
        disp("Indefinicion en alguna de las iteraciones, se retornan los valores hasta dicha iteracion");
        return
      endif
      xk=b-((b-a)/(f(b)-f(a)))*f(b);
      error=f(xk);
      
      if(abs(error)<=tol)
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
    return
  endif
  
  disp("La funcion provista no cumple el teorema de Bolzano en el intervalo especificado");
  return
end