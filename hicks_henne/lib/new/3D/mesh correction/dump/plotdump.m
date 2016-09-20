clear all;
close all;
clc;
format long;

load dumpall.txt
x = dumpall(:,2);
y = dumpall(:,3);
z = dumpall(:,4);
axis equal tight;
size(sort(unique(z)))
scatter3(x,y,z,'.');

xlabel('x\rightarrow');
ylabel('y\rightarrow');
zlabel('z\rightarrow');
