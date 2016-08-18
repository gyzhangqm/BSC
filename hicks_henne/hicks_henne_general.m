%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Hicks-Henne Bump function Parameterization                   %
%                   16 Aug 2016                                      %
%   Input: Baseline Configuration text files(Upper and lower),       %
%           Design variables                                         %
%   Output: Modified Geometry, Gradient of geometry with respect     %
%       to design variables                                          %
%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
format long;
clc;

% load airfoils/naca0012lower.txt      % Baseline configuration
% load airfoils/naca0012upper.txt
% ubase = naca0012upper(:,2);
% xu = naca0012upper(:,1);
% lbase = naca0012lower(:,2);
% xl = naca0012lower(:,1);

load airfoils/du96lower.txt      % Baseline configuration
load airfoils/du96upper.txt
ubase = du96upper(:,2);
xu = du96upper(:,1);
lbase = du96lower(:,2);
xl = du96lower(:,1);

n = 9;                      %Number of design variables
h = 1/(n+1);
bump_pos = 0+h:h:1-h;       %Bump positions
dpa = zeros(n,1);
dpb = zeros(n,1);
t_b = 4;                    %Width control parameter
for i = 1:size(ubase)
    suma = 0;
    for j = 1:n
        m = log(0.5)/log(bump_pos(j));
        suma = suma + dpa(j)*(sin(pi*xu(i)^m)^t_b);
        gradu(i,j) = sin(pi*xu(i)^m)^t_b;
    end
    unew(i) = ubase(i) + suma;
end
for i = 1:size(lbase)
    sumb = 0;
    for j = 1:n
        m = log(0.5)/log(bump_pos(j));
        sumb = sumb + dpb(j)*(sin(pi*xl(i)^m)^t_b);
        gradl(i,j) = sin(pi*xl(i)^m)^t_b;
    end
    lnew(i) = lbase(i) + sumb;
end
grad = [gradu ; gradl];

hold on;
plot(xu,unew,xl,lnew);
axis equal tight;
