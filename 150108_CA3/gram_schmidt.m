function f=gram_schmidt(A)
[row , col]=size(A);
y=zeros(row,col);
z=zeros(row,col);
y(1:row,1)=A(1:row,1)/norm(A(1:row,1));
for k=1:col-1
   z(1:row,k+1)=A(1:row,k+1);
    for i=1:k
 z(1:row,k+1)=z(1:row,k+1)-(transpose(A(1:row,k+1))*y(1:row,i))*y(1:row,i);
    end
   y(1:row,k+1)=z(1:row,k+1)/norm(z(1:row,k+1));
end
f=y;
end
