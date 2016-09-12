clear all;
close all;
clc;
format long;

load dumpalll.txt
load dumpallu.txt

xu = dumpallu(:,2);
yu = dumpallu(:,3);
zu = dumpallu(:,4);

xl = dumpalll(:,2);
yl = dumpalll(:,3);
zl = dumpalll(:,4);

size(sort(unique(zu)))
size(sort(unique(zl)))

hold on;
scatter3(xu,yu,zu,'b.');
scatter3(xl,yl,zl,'r.');
axis equal tight;

xlabel('x\rightarrow');
ylabel('y\rightarrow');
zlabel('z\rightarrow');
