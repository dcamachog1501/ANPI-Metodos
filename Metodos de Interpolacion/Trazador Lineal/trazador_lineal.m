function trazador_lineal()
 pkg load symbolic;
 clc;
 warning('off','all');
 
 points=[[1 1];[2 0.5];[4 0.25]]
 
 s=trazador_lineal_aux(points);
 
 display([s]);
 
  
end
function s= trazador_lineal_aux(points)
  syms x;
  n=length(points);
  s=sym([]);
  for(i=1:(n-1))
    x_i=points(i,1);
    y_i=points(i,2);
    x_i_1=points(i+1,1);
    y_i_1=points(i+1,2);
    s(i)=(y_i+((y_i_1-y_i)/(x_i_1-x_i))*(x-x_i));
  end
  return;
end
  