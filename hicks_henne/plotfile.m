clear all;
close all;
clc;

load updatedu.txt
load updatedl.txt
unew = updatedu(:,2);
xu = updatedu(:,1);
lnew = updatedl(:,2);
xl = updatedl(:,1);

plot(xu,unew,xl,lnew);
axis equal tight;