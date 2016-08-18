%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Hicks-Henne Polynomial Implementation       %
%       16 Aug 2016                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
format long;
clc;

N = 50;
c = 1;
r = c/2;
a = 0;
b = 2*r;
h = (b-a)/N;
ubase = zeros(N+1,1);
lbase= zeros(N+1,1);
x = a:h:b;
for i = 1:N+1
    ubase(i) = sqrt(r^2 - (x(i)-r)^2);
    lbase(i) = -sqrt(r^2 - (x(i)-r)^2);
end

dpa = zeros(5,1);
dpb = zeros(5,1);


for i = 1:N+1
    unew(i) = ubase(i) + (dpa(1)*sqrt(x(i))*(1-x(i))/exp(15*x(i))) + ...
        (dpa(2)*(sin(pi*(x(i)^0.25))^3)) + (dpa(3)*(sin(pi*(x(i)^0.757))^3))+ ...
        (dpa(4)*(sin(pi*(x(i)^1.357))^3)) + (dpa(5)*sqrt(x(i))*(1-x(i))/exp(10*x(i)));
    lnew(i) = lbase(i) + (dpb(1)*sqrt(x(i))*(1-x(i))/exp(15*x(i))) + ...
        (dpb(2)*(sin(pi*(x(i)^0.25))^3)) + (dpb(3)*(sin(pi*(x(i)^0.757))^3))+ ...
        (dpb(4)*(sin(pi*(x(i)^1.357))^3)) + (dpb(5)*sqrt(x(i))*(1-x(i))/exp(10*x(i)));
end
plot(x,unew,'*',x,lnew,'*');
hold on;
plot(x,unew,x,lnew);
axis equal tight;

