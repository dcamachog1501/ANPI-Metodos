function Descenso_Coordinado
  warning('off', 'all');
  clc;clear;
  pkg load symbolic;
  f="x*exp(-x^2-y^2)";
  variables="x y";
  iniciales="-0.5 -0.5"; %Componentes del vector inicial
  tol=1e-2;
  iterMax=2500;
  [xk,error,k]=dc(f,variables,iniciales,tol,iterMax);
  display(["El punto aproximado por el metodo es (",num2str(xk(1)),",",num2str(xk(2)) ,") con un error de ",num2str(error)," en ",num2str(k)," iteraciones"]);
end

function[xk,error,k]=dc(f,variables,iniciales,tol,iterMax)
  eval(["syms " variables ";"])
  eval(["vars=[" variables "];"]);
  eval(["xk=[" iniciales "]';"]);%Debe ser traspuestyo para poder operarlo con gk y dk
  f=sym(f);
  grdf=gradient(f,vars);
  for k=1:iterMax
    for j=1:length(vars)
      f_temp=subs(f,vars(1:end~=j), xk(1:end~=j));%(1:end~=j)->Agarrar todos los indices menos el que sea igual a j
      xk(j)=fminsearch(matlabFunction(f_temp),xk(j));
    endfor
    error=norm(double(subs(grdf,vars,xk)))
    if abs(error)<=tol
      return
    endif
  endfor
  return
end