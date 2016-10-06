clear all;
close all;
clc;
format long;

load baseline.txt
load new.txt


x = baseline(:,2);
y = baseline(:,3);
z = baseline(:,4);

xnew = new(:,2);
ynew = new(:,3);
znew = new(:,4);

figure (1);
hold on
%axis equal tight;

scatter3(x,y,z,'b.');

scatter3(xnew,ynew,znew,'r.');

xlabel('x\rightarrow');
ylabel('y\rightarrow');
zlabel('z\rightarrow');