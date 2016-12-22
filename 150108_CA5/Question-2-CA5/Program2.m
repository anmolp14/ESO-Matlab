clc
clear
syms x
fileId=fopen('CA5_q2.txt');
inputText=textscan(fileId,'%s %s','delimiter','=');
fclose(fileId);
fx=inputText{2}{1};
f=inline(fx);
x0=str2double(inputText{2}{2});
y0=str2double(inputText{2}{3});
xf=str2double(inputText{2}{4});
h=str2double(inputText{2}{5});
hmax=str2double(inputText{2}{6});
alp=str2double(inputText{2}{7});
usetol=inputText{2}{8};
gttol=inline(usetol);
tol=gttol(0);
display('1 Eulers Method')
display('2 Midpoint Method')
display('3 Runge Kutta Method')
display('4 RK45 Method')
method=input('Specify the Method Number ');
N=(xf-x0)/h;
N=N+1; 
x(1,1)=x0;
y(1,1)=y0;
use=hmax;
if method==1
    for i=2:N
       y(i,1)=(y(i-1,1)+h*f(x(i-1,1),y(i-1,1)));
       x(i,1)=(x(i-1,1)+h);
    end
    z=vpa(([x y]),5);
    z1 = double(z);
    dlmwrite('output21.txt','   x    y-euler','delimiter','')
    dlmwrite('output21.txt',z1,'-append','Delimiter','\t','precision','%.5f')
    fprintf('<<-----See "output21.txt" for answers------->>\n');
    plot(x,y)
elseif method==2
    for i=2:N
        y(i,1)=vpa(y(i-1,1)+h*f(x(i-1,1),y(i-1,1)),5);
        aver=vpa((y(i-1,1)+y(i,1))/2,5);
        x(i,1)=vpa(x(i-1,1)+h,5);
        averg=vpa((x(i-1,1)+x(i,1))/2,5);
        y(i,1)=vpa((y(i-1,1)+h*f(averg,aver)),5);
    end
    z=vpa(([x y]),5);
    plot(x,y,'-bp','DisplayName', 'Midpoint Method')
    grid on
    xlabel('x')
    ylabel('y')
    legend('-DynamicLegend','location','northwest');
    hold on; 
    z1 = double(z);
    dlmwrite('output22.txt','   x    y-midpoint','delimiter','')
    dlmwrite('output22.txt',z1,'-append','Delimiter','\t','precision','%.5f')
    fprintf('<<-----See "output22.txt" for answers------->>\n');
    
    
elseif method==3
    for i=2:N
        % All pre-defined formulas
       k1=(f(x(i-1,1),y(i-1,1)));
       k2=(f(x(i-1,1)+(h/2),y(i-1,1)+(k1*h/2)));
       k3=(f(x(i-1,1)+(h/2),y(i-1,1)+(k2*h/2)));
       k4=(f(x(i-1,1)+h,y(i-1,1)+h*k3));
       y(i,1)=(y(i-1,1)+h*(k1+2*k2+2*k3+k4)/6);
       x(i,1)=x(i-1,1)+h;
    end
    z=vpa(([x y]),5);
    z1 = double(z);
    dlmwrite('output23.txt','   x    y-RK4','delimiter','')
    dlmwrite('output23.txt',z1,'-append','Delimiter','\t','precision','%.5f')
    plot(x,y,'-go','DisplayName','RK 4th order')
    grid on
    xlabel('x')
    ylabel('y')
    legend('-DynamicLegend','location','northwest');
    fprintf('<<-----See "output23.txt" for answers------->>\n');
    
    
elseif method==4
    i=2;
    while i>=2
        j=1;
        while j==1
            k1=vpa((f(x(i-1,1),y(i-1,1))),5);
            k2=vpa((f(x(i-1,1)+(h/5),y(i-1,1)+(k1*h/5))),5);
            k3=vpa((f(x(i-1,1)+(3*h/10),(y(i-1,1)+(k1*3*h/40)+(k2*9*h/40)))),5);
            k4=vpa((f(x(i-1,1)+(3*h/5),(y(i-1,1)+(k1*3*h/10)-(k2*9*h/10)+(k3*6*h/5)))),5);
            k5=vpa((f(x(i-1,1)+(h),(y(i-1,1)-(k1*11*h/54)+(k2*h*5/2)-(k3*70*h/27)+(35*k4*h/27)))),5);
            k6=vpa((f(x(i-1,1)+(7*h/8),(y(i-1,1)+(k1*1631*h/55296)+(k2*175*h/512)+(k3*575*h/13824)+(k4*h*44275/110592)+(k5*253*h/4096)))),5);
            y5=vpa((y(i-1,1)+(2825*k1*h/27648)+(18575*k3*h/48384)+(13525*k4*h/55296)+(277*k5*h/14336)+(k6*h/4)),5);
            y4=vpa((y(i-1,1)+(37*k1*h/378)+(250*h*k3/621)+(125*k4*h/594)+(512*k6*h/1771)),5);
            e=(abs(y5-y4));
            if e<=tol
                y(i,1)=(y5);
                x(i,1)=(x(i-1,1))+(h);
                break;
            else
                kal=(h*((tol/e)^alp));
                h=kal;
            end
        end
        if x(i,1)>=xf
           break;
        end
        i=i+1;
    end
    z=vpa([x y],5);
    z1 = double(z);
    dlmwrite('output24.txt','   x    y-RK45','delimiter','')
    dlmwrite('output24.txt',z1,'-append','Delimiter','\t','precision','%.5f')
    fprintf('<<-----See "output24.txt" for answers------->>\n');
    plot(x,y)
    grid on;
    legend('RK 45');
end