clear all;
close all;
clc;
format long;

load new.txt
load baseline.txt

xnew = new(:,2);
ynew = new(:,3);
znew = new(:,4);

x = baseline(:,2);
y = baseline(:,3);
z = baseline(:,4);

figure(1)
scatter3(x,y,z,'b.');
axis equal tight;

xlabel('x\rightarrow');
ylabel('y\rightarrow');
zlabel('z\rightarrow');

figure(2)
scatter3(xnew,ynew,znew,'b.');
axis equal tight;

xlabel('x\rightarrow');
ylabel('y\rightarrow');
zlabel('z\rightarrow');