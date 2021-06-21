function lagrange()
  clc;
  warning('off','all');
  pkg load symbolic;
  
  x=[2 2.5 3 3.5 4];
  y=[4 0 -8 0 16];
  
  p=simplify(lagrangeAux(x,y))
  ezplot(p,[0 6])
end

function p=lagrangeAux(xv,y)
  syms x;
  n=length(xv);%Grado maximo del polinomio de interpolacion
  
  p=0;
  for k=1:n
    p+=y(k)*Lk(xv,n,k);
  endfor

  return
end

function L= Lk(xv,n,k)
  
  syms x;
  L=1;
  for i=1:n
    if i!=k
      L*=((x-xv(i))/(xv(k)-xv(i)));  
    endif
    
  endfor
  return;
end