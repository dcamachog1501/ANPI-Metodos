function dd_newton()
  clc;
  pkg load symbolic;
  warning('off','all');
  
  
  x=[-2 0 1];
  y=[0 1 -1];
  
  p=simplify(ddnewton_aux(x,y))
  ezplot(p,[0 6]);
  
end  

function p=ddnewton_aux(xv,y)
  syms x;
  m=length(xv);
  n=m;
  fs=[];
  for i=1:n
    f=zeros(1,m);
    for j=1:n
      if i==1
        f(j)=y(j);
        
      else
        f(j)=((fs(i-1,j+1)-fs(i-1,j))/(xv(i+j-1)-xv(j)));
      endif
    endfor
  fs=[fs;f]; 
  n-=1;  
  endfor
  
  p=0;
  
  for i=1:m
    if i==1
      p+=fs(i,1);
    else
      c=1;
      for j=1:i-1
        c*=(x-xv(j));
      endfor
      p+=fs(i,1)*c;
    endif
  endfor

  return;
 
 
end