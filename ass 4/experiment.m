
A=dlmread('input.txt',' ',0,0);
n=max(size(A));
d=input('degree=');
for i=1:d+1
    for j=1:d+1
        B(i,j)=sum(A(:,1).^((i+j)-2));
    end;
end;
B(1,1)=n;
for i=1:d+1
    C(i,1)=sum((A(:,2)).*(A(:,1).^(i-1)));
end;
X=inv(B)*C;
P=poly2sym(flip(X'));
py=vpa(P,5);
p=flip(X');
u=X';
q=polyval(p,A(:,1));
y=A(:,2);
r=1-((sum((y-q).^2))/(sum((y-mean(y)).^2)));
c=vpa(X',5);
disp(c);
disp(py);
disp(r);
plot(A(:,1),A(:,2),'square');
hold on;
h=fplot(py,[0,1]);
title('polynomial of least squares');
xlabel('x');
ylabel('y');
if (d==1)
    filename=fullfile('output.txt');
fid=fopen(filename,'wt');
fprintf(fid,'%s\n','linear');
fprintf(fid,'%s','coefficient:');
fprintf(fid,[repmat('%f\t', 1, size(u, 2)) '\n'], u');
fprintf(fid,'%s','R_sq:');
fprintf(fid,'%f\n',r);
fprintf(fid,'%s\r\n','colour:blue');
set(h,'color','b');
end;
if(d==2)
    fprintf(fid,'%s\n','quadratic');
fprintf(fid,'%s','coefficient:');
fprintf(fid,[repmat('%f\t', 1, size(u, 2)) '\n'], u');
fprintf(fid,'%s','R_sq:');
fprintf(fid,'%f\n',r);
fprintf(fid,'%s\r\n','colour:red');
    set(h,'color','r');
end;
if(d==3)
    fprintf(fid,'%s\n','cubic');
fprintf(fid,'%s','coefficient:');
fprintf(fid,[repmat('%f\t', 1, size(u, 2)) '\n'], u');
fprintf(fid,'%s','R_sq:');
fprintf(fid,'%f',r);
fprintf(fid,'%s\r\n','colour:yellow');
    set(h,'color','y');
end;
if(d==4)
    fprintf(fid,'%s\n','quaratic');
fprintf(fid,'%s','coefficient:');
fprintf(fid,[repmat('%f\t', 1, size(u, 2)) '\n'], u');
fprintf(fid,'%s','R_sq:');
fprintf(fid,'%f\n',r);
fprintf(fid,'%s','colour:green');
    set(h,'color','g');
end;
grid;