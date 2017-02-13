clear all;
close all;
clc;
format long;

load dumpbase.txt
load dumpnew.txt

n = 64;
xbase = dumpbase(:,1);
xnew = dumpnew(:,1);
ybase = dumpbase(:,2);
ynew = dumpnew(:,2);

figure(1)

scatter(xbase,ybase,'bo','DisplayName','Baseline Configuration');
%scatter(x,y,'b.');
hold on
axis equal tight;
scatter(xnew,ynew,'ro','DisplayName','Modified Configuration');
%scatter(xnew,ynew,'r.');
[h, ~] = legend('show');
xlabel('x\rightarrow');
ylabel('y\rightarrow');