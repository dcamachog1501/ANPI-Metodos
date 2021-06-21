%Metodo de Newton_Schultz para aproximar la pseudoinversa de A
function pseudoinversa
clc;
m=10;
n=5;
tol=10^-8;

%Sistema con una unica solucion
A=[1 3 4 1 2 0;
  -1 0 1 0 2 0;
   1 1 -1 -1 2 1];
b=[1 4 1]';
m=size(A)(1);
At=Newton_Schultz(A,m,tol);
xk=(At*b)';
display(["El valor aproximado para xk es ",mat2str(xk)])
end

function[xk]=Newton_Schultz(A,m,tol)
   
  alpha=max(eig(A*A'));
  xk=A'/alpha^2;% Valor incial de X.
  I=eye(m);
  for k=1:15
    xk=xk*(2*I-A*xk);
    
    frobeniusmat=A*xk*A-A;
    error=norm(frobeniusmat,'fro');
    if error<=tol
      break;
    endif;
  end;
  return 
  
end