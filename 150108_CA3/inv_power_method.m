function [ Y, Mf, totitr ] = inv_power_method (len , A, maxitr, relerr )
   A=inv(A);
    for j=1:len
       guess(j,1)=1; 
    end
    Mi=0;
    for i=1:maxitr
        Y=A*guess;
        Mf=max(abs(Y));
        if(abs(((Mi-Mf)/Mf))*100<relerr||i==maxitr)
            totitr=i;
            break
        end
        guess=Y/Mf;
        Mi=Mf;
    end
    Mf= 1.0/Mf;
    Y = Y / norm(Y);
    fileId = fopen('output inverse power.txt','w');
    fprintf(fileId,'Eigenvalue \n %.4f \n \n', Mf );
    fprintf(fileId,'Eigenvector \n ' );
    fprintf(fileId,' \n %.4f ', Y );
    fprintf(fileId,'\n \n Number of iterations \n %d', totitr );
    disp('eigen value');
    disp(Mf);
    disp('eigen vector');
    disp(Y);
    disp('number of iterations');
    disp(totitr);
    disp('the output is also displayed in output inverse power.txt file');