function test
  tridiagonal(1,2,3,4)
end
function[a_mat]=tridiagonal(a,b,c,n)
  clc;
  a_mat=zeros(n,n);
  for i=1:n
    if i==1
      a_mat(i,i+1)=b;
    elseif i==n-1
      a_mat(i,i-1)=c;
    else
      a_mat(i,i+1)=b;
      a_mat(i,i-1)=c;
    endif
    a_mat(i,i)=a;
  endfor
  return
end