clear all;
close all;
clc;
format long;

load baseline.txt
load new.txt


x = baseline(:,2);
y = baseline(:,3);
z = baseline(:,4);

xnew = new(:,2);
ynew = new(:,3);
znew = new(:,4);

figure (1);
hold on
scatter3(x,y,z,'b.','DisplayName','Baseline Configuration');
scatter3(xnew,ynew,znew,'r.','DisplayName','Optimized Configuration');
xlabel('x\rightarrow');
ylabel('y\rightarrow');
zlabel('z\rightarrow');
[h, ~] = legend('show');



% j = 1;
% for i = 1:size(x)
%     if (y(i) <= ypanel+eps) && (y(i) >= ypanel-eps)
%         oldx(j) = x(i);
%         oldy(j) = y(i);
%         oldz(j) = z(i);
%         upx(j) = xnew(i);
%         upy(j) = ynew(i);
%         upz(j) = znew(i);
%         j = j + 1;
%     end
% end
% 
% hale = delaunayTriangulation(oldx',oldy',oldz');
% hosa = delaunayTriangulation(upx',upy',upz');
% %faceColor  = [0.6875 0.8750 0.8984];
% faceColor = [0 0 1];
% faceColor2 = [1 0 0];
% 
% hold on;
% tetramesh(hale,'FaceColor',faceColor,'FaceAlpha',0.3,'edgecolor','none');
% tetramesh(hosa,'FaceColor',faceColor2,'FaceAlpha',0.3,'edgecolor','none');
% axis equal tight;