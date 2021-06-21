function Gradiente_Conjugado_No_Lineal
  clc;clear;
  pkg load symbolic;
  f="(x+9)^2+y^2+5";
  variables="x y";
  iniciales="3 0"; %Componentes del vector inicial
  tol=1e-5;
  iterMax=2500;
  [xk,error,k]=gcnl(f,variables,iniciales,tol,iterMax);
  display(["El punto aproximado por el metodo es (",num2str(xk(1)),",",num2str(xk(2)) ,") con un error de ",num2str(error)," en ",num2str(k)," iteraciones"]);
end

function[xk,error,k]=gcnl(f,variables,iniciales,tol,iterMax)
  warning('off', 'all');
  eval(["syms " variables ";"])
  eval(["vars=[" variables "];"]);
  eval(["xk=[" iniciales "]';"]);%Debe ser traspuesto para poder operarlo con gk y dk
  f=sym(f);
  grdf=gradient(f,vars);
  grs=sym(grdf);
  gk=double(subs(grs,vars,xk));
  dk=-gk;
  k=0;error=0;
  
  while k<iterMax
    delta=rand(1);
    ak=0.5;
    aux=xk+ak*dk;
    izq=double(subs(f,vars,aux));
    der=delta*ak*gk'*dk;
    while izq>der
      ak=ak/2;
      aux=xk+ak*dk;
      izq=double(subs(f,vars,aux)-subs(f,vars,xk));
      der=delta*ak*gk'*dk;
    endwhile
    xk=xk+ak*dk;
    error=norm(double(subs(grs,vars,xk)))
    if abs(error)<=tol
      error=abs(error);
      return
    endif
    gkp1=double(subs(grs,vars,xk));
    beta=norm(gkp1)^2/norm(gk)^2;
    gk=gkp1;
    dk=-gk+beta*dk;
    k=k+1;
  endwhile
  error=abs(error);
  return;
end