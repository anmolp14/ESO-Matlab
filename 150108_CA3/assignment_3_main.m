
fileID = fopen('input.txt','r');
n = textscan(fileID,'%d',1,'Delimiter','\n');
M = dlmread('input.txt',' ',[1 0 n{1} n{1}-1 ]);
max_iter = dlmread ('input.txt',' ',[n{1}+1 0 n{1}+1 0]);
max_rel_error = dlmread ('input.txt',' ',[n{1}+2 0 n{1}+2 0]) ;
shift_factor = dlmread ('input.txt',' ',[n{1}+3 0 n{1}+3 0]) ;
no = input('enter the method number you want to use \n 1.) Direct power method \n 2.) Inverse power method \n 3.) Shifted power method \n 4.) QR method \n');
 
if no ==1
    direct_power_method(n{1},M,max_iter,max_rel_error);
   
elseif no ==2
    inv_power_method(n{1},M,max_iter,max_rel_error);
    
elseif no==3    
    Shifted_power(n{1},M,max_iter,max_rel_error,shift_factor);
    
elseif no==4
    QR_method(n{1},M, max_iter,max_rel_error);
else 
    disp('Error:please enter a number between 1-4');
end    



