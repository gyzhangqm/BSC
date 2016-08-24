clear all;
close all;
clc;
format long;

load updatedu.txt
load updatedl.txt
Nu = 72;
Nl = 72;
panels = 5;


xu = updatedu(:,1);
xl = updatedl(:,1);
yu = updatedu(:,2);
yl = updatedl(:,2);
zu = updatedu(:,3);
zl = updatedl(:,3);

U = delaunayTriangulation(xu,yu,zu);
L = delaunayTriangulation(xl,yl,zl);

faceColor  = [0.6875 0.8750 0.8984];
figure
hold on;
tetramesh(L,'FaceColor',faceColor,'FaceAlpha',0.3,'edgecolor','black');
tetramesh(U,'FaceColor',faceColor,'FaceAlpha',0.3,'edgecolor','black');

for i = 1:Nu*panels
    scatter3(updatedu(i,1), updatedu(i,2), updatedu(i,3),'r.');
    scatter3(updatedl(i,1), updatedl(i,2), updatedl(i,3),'r.');
end
axis equal tight;

axis equal tight;