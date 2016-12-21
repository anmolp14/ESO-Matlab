display=('enter the method name:- \n L for linear spline \n Q for quadratic spline \n C for cubic spline \n N for knot-a-knot spline \n P for periodic splne \n CL for clamped spline');
method=input(display,'s');
disp(method);
x=dlmread('inputx.txt',' ',0,0);
y=dlmread('inputy.txt',' ',0,0);
xc=dlmread('textc.txt',' ',0,0);
n=max(size(x));
if(method=='L')
    for i=1:(max(size(x))-1);
    q(i)=((xc(i)-x(i+1))/(x(i)-x(i+1)))*y(i)+((xc(i)-x(i))/(x(i+1)-x(i)))*y(i+1);
end;
d=xc';
e=q';
R=[d,e];
disp('colour:red');
for i=1:(max(size(x))-1);
    syms z;
    p(i)=((z-x(i+1))/(x(i)-x(i+1)))*y(i)+((z-x(i))/(x(i+1)-x(i)))*y(i+1);
end;
for j=1:(max(size(x))-1);
     ylim([-1,1.5]);
    fplot(p(j),[x(j),x(j+1)],'r');
    hold on;
    title('fitted splines polynomial graph');
end;
end;
if(method=='Q');
    c(2)=0;
        A=[1 x(1);1 x(2)];
        d=[y(1);y(2)];
        e=inv(A)*d;
        for i=1:n
            a(2)=e(1);
            b(2)=e(2);
        end;
for i=3:n
        B=[0 1 2*x(i-1);1 x(i-1) (x(i-1))^2;1 x(i) (x(i))^2];
        f=[b(i-1)+2*(x(i-1)*c(i-1));y(i-1);y(i)];
        g=inv(B)*f;
        a(i)=g(1);
        b(i)=g(2);
        c(i)=g(3);
end;
for i=1:n-1
    q(i)=a(i+1)+b(i+1)*xc(i)+c(i+1)*xc(i).^2;
end;
e=q';
d=xc';
R=[d,e];
disp(R);
disp('colour:blue');
for i=1:n-1
    syms z;
    p(i)=a(i+1)+b(i+1)*z+c(i+1)*z^2;
end;
for j=1:(max(size(x))-1)
    xlabel('x');
    ylabel('y');
    ylim([-1,1.5]);
    fplot(p(j),[x(j),x(j+1)],'b');
    hold on;
end;
end;
if(method=='C')
    for i=1:n
    if(i==1)
        h(i)=0;
        end;
    if(i~=1)
    h(i)=(x(i)-x(i-1));
    end;
end;
for i=2:n
    g(i)=((y(i)-y(i-1))/h(i));
end;
T=zeros(n-2,n-2);
for i=1:n-3
    T(i,i+1)=h(i+2);
    T(i+1,i)=h(i+2);
end
for i=1:n-2
    T(i,i)=2*(h(i+1)+h(i+2));
end
for j=1:n-2
    b(j)=6*(g(j+2)-g(j+1));
end
Y=inv(T)*b';
for i=1:n
    s(i)=0;
end
for i=1:n-2
    s(i+1)=Y(i);
end
for i=2:n
  A(i)=(s(i)/(6*h(i)));
  B(i)=(s(i-1)/(6*h(i)));
  C(i)=((y(i)/h(i))-((s(i)*h(i))/6));
  D(i)=((y(i-1)/h(i))-((s(i-1)*h(i))/6));  
end
for i=1:(max(size(x))-1)
q(i)=A(i+1)*((xc(i)-x(i))^3)-B(i+1)*((xc(i)-x(i+1))^3)+C(i+1)*(xc(i)-x(i))-D(i+1)*(xc(i)-x(i+1));
end;
e=q';
d=xc';
R=[d,e];
disp(R);
disp('colour:green');
for i=1:(max(size(x))-1)
    syms y;
p(i)=A(i+1)*((y-x(i))^3)-B(i+1)*((y-x(i+1))^3)+C(i+1)*(y-x(i))-D(i+1)*(y-x(i+1));
pj(i)=vpa(p(i),4);
end;
for j=1:(max(size(x))-1)
    ylim([-1,1.5]);
    fplot(pj(j),[x(j),x(j+1)],'g');
    hold on;
end;
end;
if(method=='P')
    for i=1:n
    if(i==1)
        h(i)=0;
        end;
    if(i~=1)
    h(i)=(x(i)-x(i-1));
    end;
end;
for i=2:n
    g(i)=((y(i)-y(i-1))/h(i));
end;
T=zeros(n,n);
T(1,1)=(5*h(2))/6;
T(1,2)=h(2)/6;
T(1,n-1)=h(n)/6;
T(1,n)=h(n)/3;
for j=2:n-1
    T(j,j+1)=h(j+1);
    T(j,j-1)=h(j);
    T(j,j)=2*(h(j)+h(j+1));
end;
T(n,1)=1;
T(n,n)=-1;
b(1)=((y(2)-y(n))/h(2))+((y(n-1)-y(n))/h(n));
b(n)=0;
for j=2:n-1
    b(j)=6.*(g(j+1)-g(j));
end;
Y=inv(T)*b';
for i=1:n
    s(i)=Y(i);
end;
for i=2:n
  A(i)=(s(i)/(6*h(i)));
  B(i)=(s(i-1)/(6*h(i)));
  C(i)=((y(i)/h(i))-((s(i)*h(i))/6));
  D(i)=((y(i-1)/h(i))-((s(i-1)*h(i))/6));  
end
for i=1:(max(size(x))-1)
q(i)=A(i+1)*((xc(i)-x(i))^3)-B(i+1)*((xc(i)-x(i+1))^3)+C(i+1)*(xc(i)-x(i))-D(i+1)*(xc(i)-x(i+1));
end;
e=q';
d=xc';
R=[d,e];
disp(R);
disp('colour:green--');
for i=1:(max(size(x))-1)
    syms z;
p(i)=A(i+1)*((z-x(i))^3)-B(i+1)*((z-x(i+1))^3)+C(i+1)*(z-x(i))-D(i+1)*(z-x(i+1));
pj(i)=vpa(p(i),4);
end;
for j=1:(max(size(x))-1)
    fplot(pj(j),[x(j),x(j+1)],'g--');
    hold on;
end;
end;
if(method=='CL')
    a=input('si=');
c=input('sf=');
for i=1:n
    if(i==1)
        h(i)=0;
        end;
    if(i~=1)
    h(i)=(x(i)-x(i-1));
    end;
end;
for i=2:n
    g(i)=((y(i)-y(i-1))/h(i));
end;
T=zeros(n,n);
T(1,1)=(5*h(2))/6;
T(1,2)=h(2)/6;
for j=2:n-1
    T(j,j+1)=h(j+1);
    T(j,j-1)=h(j);
    T(j,j)=2*(h(j)+h(j+1));
end;
T(n,n-1)=h(n)/6;
T(n,n)=h(n)/3;
b(1)=((y(2)-y(1)-(h(2)*a)))/h(2);
b(n)=(c+(y(n-1)/h(n-1))-(y(n)/h(n)));
for j=2:n-1
    b(j)=6*(g(j+1)-g(j));
end;
Y=inv(T)*b';
for i=1:n
    s(i)=Y(i);
end;
for i=2:n
  A(i)=(s(i)/(6*h(i)));
  B(i)=(s(i-1)/(6*h(i)));
  C(i)=((y(i)/h(i))-((s(i)*h(i))/6));
  D(i)=((y(i-1)/h(i))-((s(i-1)*h(i))/6));  
end
for i=1:(max(size(x))-1)
q(i)=A(i+1)*((xc(i)-x(i))^3)-B(i+1)*((xc(i)-x(i+1))^3)+C(i+1)*(xc(i)-x(i))-D(i+1)*(xc(i)-x(i+1));
end;
e=q';
d=xc';
R=[d,e];
disp(R);
disp('colour:red--');
for i=1:(max(size(x))-1)
    syms y;
p(i)=A(i+1)*((y-x(i))^3)-B(i+1)*((y-x(i+1))^3)+C(i+1)*(y-x(i))-D(i+1)*(y-x(i+1));
pj(i)=vpa(p(i),4);
end;
for j=1:(max(size(x))-1)
    ylim([-1,1.5]);
    fplot(pj(j),[x(j),x(j+1)],'r--');
    hold on;
end;
end;
if(method=='N')
    
for i=1:n
    if(i==1)
        h(i)=0;
        end;
    if(i~=1)
    h(i)=(x(i)-x(i-1));
    end;
end;
for i=2:n
    g(i)=((y(i)-y(i-1))/h(i));
end;
T=zeros(n,n);
T(1,1)=h(3);
T(1,2)=-(h(2)+h(3));
T(1,3)=h(2);
T(n,n-2)=h(n);
T(n,n-1)=-(h(n-1)+h(n));
T(n,n)=h(n-1);
for j=2:n-1
    T(j,j+1)=h(j+1);
    T(j,j-1)=h(j);
    T(j,j)=2*(h(j)+h(j+1));
end;
b(1)=0;
b(n)=0;
for j=2:n-1
    b(j)=6.*(g(j+1)-g(j));
end;
Y=inv(T)*b';
for i=1:n
    s(i)=Y(i);
end;
for i=2:n
  A(i)=(s(i)/(6*h(i)));
  B(i)=(s(i-1)/(6*h(i)));
  C(i)=((y(i)/h(i))-((s(i)*h(i))/6));
  D(i)=((y(i-1)/h(i))-((s(i-1)*h(i))/6));  
end
for i=1:(max(size(x))-1)
q(i)=A(i+1)*((xc(i)-x(i))^3)-B(i+1)*((xc(i)-x(i+1))^3)+C(i+1)*(xc(i)-x(i))-D(i+1)*(xc(i)-x(i+1));
end;
e=q';
d=xc';
R=[d,e];
disp(R);
disp('colour:b*');
for i=1:(max(size(x))-1)
    syms y;
p(i)=A(i+1)*((y-x(i))^3)-B(i+1)*((y-x(i+1))^3)+C(i+1)*(y-x(i))-D(i+1)*(y-x(i+1));
pj(i)=vpa(p(i),4);
end;
for j=1:(max(size(x))-1)
    ylim([-1 1.5]);
    fplot(pj(j),[x(j),x(j+1)],'b*');
    hold on;
end;
end;


    