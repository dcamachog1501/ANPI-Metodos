function runge_kutta_4
  pkg load symbolic;
  warning('off','all');
  clc;
  
  
  f='y-x^2+1';
  y0=0.5;
  a=0;b=2;
  n=11;
  [xv,yv]=RK4(f,a,b,y0,n);
  
  %Graficacion del polinomio de interpolacion y los puntos (xv,yv)
  stem(xv,yv,'LineStyle','none');


end

function [xv,yv]=RK4(f,a,b,y0,n)
  
  h=(b-a)/(n-1);
  xv=a:h:b;
  f=matlabFunction(sym(f));
  yv=zeros(1,n);
  yv(1)=y0;
  for k=1:n-1
    k1=f(xv(k),yv(k));
    k2=f(xv(k)+h/2,yv(k)+((h*k1)/2));
    k3=f(xv(k)+h,yv(k)+h*(2*k2-k1));
    k4=f(xv(k)+h,yv(k)+h*k3);
    yv(k+1)=yv(k)+(k1+2*k2+2*k3+k4)*h/6;
  endfor
  return
  
end