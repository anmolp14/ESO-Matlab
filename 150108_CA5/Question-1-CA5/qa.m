function [i n x] = qa( str, a, b, x, tol, n)

f = str2func(str);
c = (a+b)/2;
h=(b-a)/2;
ih = h*(f(a)+4*f(c)+f(b))/3;
ihb2 = h/6*(f(a)+4*f((a+c)/2)+f(c)) + h/6*(f(c)+4*f((c+b)/2)+f(b));
e = abs(ihb2-ih);
x(n) = (a+c)/2;
x(n+1) = (c+b)/2;
if(e>=tol)
    [ih1,n,x] = qa(str,a,c,x,tol,n+2);
    [ih2,n,x] = qa(str,c,b,x,tol,n+2);
    i = ih1+ih2;
else
    i = 16/15*ihb2 - 1/15*ih;
end
end



