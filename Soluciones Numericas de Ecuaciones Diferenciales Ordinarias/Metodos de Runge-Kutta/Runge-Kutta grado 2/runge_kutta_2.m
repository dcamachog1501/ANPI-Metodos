function runge_kutta_2
  pkg load symbolic;
  warning('off','all');
  clc;
  
  
  f='y-x^2+1';
  y0=0.5;
  a=0;b=2;
  n=11;
  [xv,yv]=RK2(f,a,b,y0,n);
  
  %Graficacion del polinomio de interpolacion y los puntos (xv,yv)
  stem(xv,yv,'LineStyle','none');


end

function [xv,yv]=RK2(f,a,b,y0,n)
  
  h=(b-a)/(n-1);
  xv=a:h:b;
  f=matlabFunction(sym(f));
  yv=zeros(1,n);
  yv(1)=y0;
  for k=1:n-1
    k1=f(xv(k),yv(k));
    k2=f(xv(k)+h/2,yv(k)+((h*k1)/2));
    yv(k+1)=yv(k)+h*k2;
  endfor
  return
  
end