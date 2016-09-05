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
hold on;
axis equal tight;
for i = 1:n
    scatter(xbase(i),ybase(i));
end
figure(2)
hold on;
axis equal tight;
for i = 1:n
    scatter(xnew(i),ynew(i));
end