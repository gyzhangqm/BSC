clear all;
close all;
clc;
format long;

load new.txt
load baseline.txt

xnew = new(:,2);
ynew = new(:,3);


x = baseline(:,2);
y = baseline(:,3);


figure(1)
scatter(x,y,'b.');
axis equal tight;

xlabel('x\rightarrow');
ylabel('y\rightarrow');
zlabel('z\rightarrow');

figure(2)
scatter(xnew,ynew,'b.');
axis equal tight;

xlabel('x\rightarrow');
ylabel('y\rightarrow');
zlabel('z\rightarrow');