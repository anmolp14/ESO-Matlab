clc
clear
syms x
fileId=fopen('CA5_q1.txt');
inputText=textscan(fileId,'%s %s','delimiter','=');
fclose(fileId);
st=inputText{2}{1};
str = strcat('@(x)',st);
f = str2func(str);
a=str2double(inputText{2}{2});
b=str2double(inputText{2}{3});
tl=(inputText{2}{4});
tol=inline(tl);
tol=tol(0);
x(1)=a;
x(2)=b;
x(3)=(a+b)/2;
[i,n,x] = qa(str,a,b,x,tol,4);
x = sort(x);
for j=1:n+1
    y(j)=f(x(j));
end;
plot(x, y, '-o');
xlabel('x');
st = strcat('f(x)=',st);
ylabel(st);
title('Rel. Approx. Error vs. Iteration Number');
grid;
file = fopen('output1.txt','w');
fprintf(file,'I = %.4f\n', i);
fprintf(file,'n = %d',n);
fprintf('<<-----See "output1.txt" for answers------->>\n');
fclose(file);


    
 

