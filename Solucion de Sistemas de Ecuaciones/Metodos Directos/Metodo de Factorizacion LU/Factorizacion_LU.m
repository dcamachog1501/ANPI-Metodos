
function Factorizacion_LU
  clc;
  a=[4,-2,1;20,-7,12;-8,13,17];
  b=[11;70;17];
  [xk]=fact_lu(a,b);
  if xk!=inf
    xk_str=[];
      for i=1:size(xk,2)  
        xk_str=[xk_str "x" num2str(i) "=" num2str(xk(i)) " "];
      endfor
      disp(["Las soluciones del sistema de ecuaciones son: "]);
      disp(xk_str);
  endif
end
function[xk]=fact_lu(a,b)
  [l,u]=lu(a);
  if l!=inf
    [yk]=sust_adelante(l,b);
    [xk]=sust_atras(u,yk');%yk debe ser traspuesto para que sea un vector columna.
  else
    xk=inf;
  endif
  return
  
end
function[l,u]=lu(a)
  [n,m]=size(a);
  if m==n
    for i=1:n%Verificacion de existencia y unicidad
      temp=a(1:i,1:i);
      if det(temp)==0
        l=inf;
        u=0;
        disp("Una de las submatrices de a no es invertible");
        return
      endif
    endfor
    
    l=eye(n);
    for k=1:n
      for i=k+1:n
        l(i,k)=a(i,k)/a(k,k);
        a(i,k)=0;
        for j=k+1:n
          a(i,j)=a(i,j)-l(i,k)*a(k,j);%Transformando a en triangular superior
        endfor
      endfor
    endfor
    u=a;
    return  
  else
    l=inf;
    u=0;
    disp("La matriz provista no es cuadrada");
    return
  endif
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
function[xk]=sust_adelante(a,b)
  n=size(b,1);%numero de filas en la matriz
  xk=zeros(1,n);%array de 1xn dimensiones lleno de n 0's
  for i=1:n
    suma=0;
    for j=1:i
      suma=suma+a(i,j)*xk(j);
    endfor
    xi=(1/a(i,i))*(b(i)-suma);
    xk(i)=xi;
  endfor
  return
end