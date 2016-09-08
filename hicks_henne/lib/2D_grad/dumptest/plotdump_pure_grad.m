clear all;
close all;
clc;
format long;

load dump.txt


n = 64;
xbase = dump(:,1);
ybase = dump(:,2);

figure(1)
hold on;
axis equal tight;
for i = 1:n
    scatter(xbase(i),ybase(i),'r');
end
