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



scatter(x,y,'b.','DisplayName','Baseline Configuration');
%scatter(x,y,'b.');
hold on
%axis equal tight;
scatter(xnew,ynew,'r.','DisplayName','Optimized Configuration');
%scatter(xnew,ynew,'r.');
[h, ~] = legend('show');