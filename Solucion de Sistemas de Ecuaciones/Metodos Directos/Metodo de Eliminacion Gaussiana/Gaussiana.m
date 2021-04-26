function Gaussiana
  clc;
  warning('off', 'all');
  [xk]=gaussiana([2,-6,12,16;1,-2,6,6;-1,3,-3,-7;0,4,3,-6],[70;26;-30;-26]);
  if xk(1)!=inf
    xk_str=[];
    for i=1:size(xk,2)  
      xk_str=[xk_str "x" num2str(i) "=" num2str(xk(i)) " "];
    endfor
    disp(["Las soluciones del sistema de ecuaciones son: "]);
    disp(xk_str);
  else
    disp("El determinante de la matriz provista es 0, por lo cual esta no posee inversa, lo cual no permite utilizar este metodo");
  endif
end
function[xk]=gaussiana(a,b)
  n=size(b,1);
  xk=zeros(1,n);
  if det(a)!=0
    [at,bt]=triang_sup(a,b);
    [xk]=sust_atras(at,bt);
  else
    xk(1)=inf;
  endif
  return
end
function[xk]=sust_atras(a,b)
  n=size(b,1);%numero de filas en la matriz
  xk=zeros(1,n);%array de 1xn dimensiones lleno de n 0's
  for i=n:-1:1
    suma=0;
    for j=i+1:n
      suma=suma+a(i,j)*xk(j);
    endfor
    xi=(1/a(i,i))*(b(i)-suma);
    xk(i)=xi;
   
  endfor
  return
end
function [at,bt]=triang_sup(a,b)
  %b debe ser un vector columna para hacer la matriz aumentada 
  n=size(b,1);
  a_aumt=[a b];
  for k=1:n-1
    for i=k+1:n
      mik=a_aumt(i,k)/a_aumt(k,k);
      for j=k:n+1
        a_aumt(i,j)=a_aumt(i,j)-mik*a_aumt(k,j);
      endfor
    endfor
  endfor
  at=a_aumt(:,1:n);
  bt=a_aumt(:,n+1);
  return
  
  
end