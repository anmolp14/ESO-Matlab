function [phi , iter  ] = QR_method( len , A, maxitr, relerr )
x=zeros(len,1);
    phi=zeros(len,1);
    iter=1;
    rel= 100;
    
    %temporary A whose value changes in loop
    tempA=A;
    while((iter<=maxitr && rel>=relerr)||iter==1)
    
    %orthogonalize A
    Q=gram_schmidt(tempA);
    
    %Q inverse = Q transpose ; R = Q inverse * A
    R=transpose(Q)*tempA;
    
    % getting the diagonal elements
    x=diag(R*Q-tempA)./diag(R*Q);
    x=abs(x);
    
    % getting relative error
    rel=max(x)*100;
    
    % update value to R*Q
    tempA=R*Q; 
    iter=iter +1;
    end
    
    % remove extra addition to iter value
    iter=iter-1;
    phi=diag(tempA);
    fid = fopen('output QR method.txt', 'wt');
    fprintf(fid,'Eigenvalues \n\n');
    fprintf(fid,'%.4f \n',phi);
    fprintf(fid,'\n\nIterations \n\n');
    fprintf(fid,'%d',iter);
    disp('output is displayed in the output QR method.txt file');
end
        



