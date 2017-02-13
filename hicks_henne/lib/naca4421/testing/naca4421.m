clear;
close;
clc;
format long;
A = load('test.dat');
%scatter(A(:,2),A(:,3));

B = load('naca.fix.nod');
for i = 1:size(B,1)
    %if (B(i,2) == 3) 
        scatter(A(B(i,1),2),A(B(i,1),3));
        hold on;
    %end
end
axis equal tight;
