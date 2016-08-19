clear all;
close all;
clc;
format long;

load updatedu.txt
load updatedl.txt
Nu = 72;
Nl = 72;
panels = 3;
figure(1);
hold on;
for i = 1:Nu*panels
    scatter3(updatedu(i,1), updatedu(i,2), updatedu(i,3),'r.');
    scatter3(updatedl(i,1), updatedl(i,2), updatedl(i,3),'r.');
end
axis equal tight;
