
function bfgs
  pkg load symbolic;
  warning('off','all');
  clc;
  f=5;
  
  %Definicion al azar del vector xk inicial.
  x0=rand(1,f)
  
  B0=eye(f);
  [xk,error,k]=bfgs_aux(f,x0,B0,2500,1e-5);
  xk=xk';
  values="(";
  for i=1:size(xk)(2)+1
    if i==size(xk)(2)
      values=[values,num2str(xk(i)),")"];
      break;
    endif
    values=[values,num2str(xk(i)),","];
  endfor
  display(["El punto aproximado por el metodo es ",values," con un error de ",num2str(error)," en ",num2str(k)," iteraciones"]);
  
 
% Funcion auxiliar que implementa la funcionalidad del algoritmo BFGS
% Entradas : 
%      f : Numero de variables para la funcion a generar (Explicado más adelante)
%      x0 : Conjunto de valores xk iniciales que se utilizaran en la primera iteracion.
%      B0: matriz simetrica definida positiva de nxn inicial.
%      iterMax : numero maximo de iteraciones
%      tol: Tolerancia de error maxima para la aproximacion.
% Salida : 
%      xk : Vector con las coordenadas del valor minimo en la funcion.
%      error: Error de la aproximacion
%      k: Numero de iteraciones
function[xk,error,k]=bfgs_aux(f,x0,B0,iterMax,tol)

  %Funcion escogida: Sumatoria de xi^10 con i de 1 a n.
  %Por la naturaleza de la funcion, se decidio que en ves de pasar la funcion como parametro se pasaria el numero n de variables
  %por lo cual el primer paso del algoritmo es generar la funcion de de la sumatoria de n variables elevadas a la 10 con el siguiente bloque de codigo.
  
  faux="";
  vars='';
  for i=1:f+1
    if i==f
      faux=[faux,"x",num2str(i),"^10"];
      vars=[vars,"x",num2str(i)];
      break;
    endif
    faux=[faux,"x",num2str(i),"^10","+"];
    vars=[vars,"x",num2str(i)," "];
  endfor
  eval(["syms " vars ";"]);
  eval(["vars=[" vars "];"]);
  xk=x0';
  f=sym(faux);
  
  % Inicio del bloque propio del algoritmo BFGS
  
  %Definicion del gradiente de la funcion generada.
  grdf=gradient(f,vars);
  
  %Representacion simbolica del gradiente.
  g=sym(grdf);
  
  %Variable que almacena la matriz Bk
  Bk=B0;
  
  %Variable que almacena el error del resultado aproximado.
  error=inf;
  e=[];
  
  %Valores sugeridos para las constantes omega1,omega2,alpha y epsilon
  omega1=1e-4;
  omega2=0.9;
  alpha=5;
  epsilon=9;
 
 %Ciclo for que realiza el metiodo iterativo
  for k=0:iterMax
    
    %Definicion de la direccion BFGS pk
    pk=double(-inv(Bk)*subs(g,vars,xk));
    
    %Definicion de los valores necesarios para encontrar el lambda ideal para la iteracion.
    lambda=1;
    lambdamin=0;
    lambdamax=inf;
    
    %Valores de g y f para la iteracion k en la que se encuentra el ciclo
    gk=double(subs(g,vars,xk));
    fk=double(subs(f,vars,xk));
    
    %Ciclo while que permite encontrar el valor ideal de lamda para la iteracion utilizando una variacion del metodo Wolfe
    while true
      if double(subs(f,vars,xk+lambda*pk))>fk+omega1*lambda*gk'*pk
        lambdamax=lambda;
        lambda=(lambdamax+lambdamin)/2;
      elseif double(subs(g,vars,xk+lambda*pk)'*pk)<omega2*gk'*pk
        lambdamin=lambda;
        if lambdamax==inf
          lambda=lambda*2;
        else
          lambda=(lambdamax+lambdamin)/2;
        endif
      else
        break;
      endif
    endwhile
    
    %Definicion del xk de la siguiente iteracion.
    xk_1=double(xk+lambda*pk);
    
    %Definicion de los valores necesarios para encontrar la matriz Bk de la siguiente iteracion.
    sk=double(xk_1-xk);
    yk=double(subs(g,vars,xk_1)-gk);
    
    %Bloque condicional para determinar cual valor asignar al siguiente Bk.
    if ((yk'*sk)/(norm(sk))^2)>=(epsilon*(norm(gk))^alpha)
      Bk=double(Bk-(Bk*sk*sk'*Bk)/(sk'*Bk*sk)+(yk*yk')/(yk'*sk));
    endif
    
    %Calculo del error para esta iteracion
    error=norm(double(subs(g,vars,xk_1)))
    
    %Condicion de parada
    if error<=tol
      e(end+1)=error;
      break;
    endif
    xk=xk_1;
    e(end+1)=error;
  endfor
  
  plot(1:k+1,e,"--o")
  xlabel('Cantidad Iteraciones')
  ylabel('Error |f(Xk)|')
  title('Gráfico Iteraciones vs Error')
  
  return;
