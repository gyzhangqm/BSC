clear all;
close all;
clc;
format longe;
load plotdata.txt
x = plotdata(:,1);
l = plotdata(:,2);
d = plotdata(:,3);
a = plotdata(:,4);

% plot(x,l,'r');
% hold on;
% plot(x,d,'b');
plot(x,d)

