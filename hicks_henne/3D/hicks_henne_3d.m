%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Hicks-Henne Bump function Parameterization                   %
%                   16 Aug 2016                                      %
%             3D Wing parameterization                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
format long;

% Baseline Config Parameters Onera M6
panels = 5;
hh = 5;
t_b = 4;
span = 1196.3;
tspan = 2734.5;
xspan = 1578.7;
sweep = 26.7*pi/180;
leangle = 60*pi/180;
teangle = 74.219120*pi/180;
% Baseline Config Parameters Onera M6

load airfoils/oneraup.txt
load airfoils/oneradown.txt
load dp/dispx.txt
load dp/dispy.txt
load dp/dp.txt
load dp/scale.txt
load dp/twist.txt
ubase = oneraup(:,2);
xu = oneraup(:,1);
lbase = oneradown(:,2);
xl = oneradown(:,1);


% Var Initialisations
%ypanelpos = zeros(panels,1);
xpanelpos = zeros(panels,1);
chordlen = zeros(panels,1);
ubasepanel = zeros(panels,size(ubase,1));
lbasepanel = zeros(panels,size(lbase,1));
unewpanel = zeros(panels,size(ubase,1));
lnewpanel = zeros(panels,size(lbase,1));
xubasepanel = zeros(panels,size(ubase,1));
xlbasepanel = zeros(panels,size(lbase,1));
xunewpanel = zeros(panels,size(ubase,1));
xlnewpanel = zeros(panels,size(lbase,1));
dpa = zeros(panels,hh);
dpb = zeros(panels,hh);
twist = zeros(panels,1);
ypanelpos = 0:span/(panels-1):span;
h = 1/(hh+1);
bump_pos = 0+h:h:1-h;
% Var Initialisations


% Reading Design parameters from dp folder
dispx = dispx(:,1);
dispy = dispy(:,1);
scale = scale(:,1);
twistnew = twist(:,1);
dpas = dp(:,1);
dpbs = dp(:,2);
for i = 1:panels
    for j = 1:hh
        dpa(i,j) = dpas((i-1)*hh + j);
        dpb(i,j) = dpbs((i-1)*hh + j);
    end
end
% Reading Design parameters from dp folder


for i = 1:panels
    % Baseline Config Generation start
    xpanelpos(i) = (tspan-ypanelpos(i))/tan(0.5*pi-sweep);
    chordlen(i) = ((tspan-ypanelpos(i))/tan(leangle)) - ((tspan-ypanelpos(i))/tan(teangle));
        %Scaling
    ubasepanel(i,:) = chordlen(i).*ubase;
    lbasepanel(i,:) = chordlen(i).*lbase;
    xubasepanel(i,:) = chordlen(i).*xu;
    xlbasepanel(i,:) = chordlen(i).*xl;
        %Translating
    xubasepanel(i,:) = xubasepanel(i,:) - (tspan-ypanelpos(i))/tan(leangle);
    xlbasepanel(i,:) = xlbasepanel(i,:) - (tspan-ypanelpos(i))/tan(leangle);
        % Twisting
    xubasepanel(i,:) = (xubasepanel(i,:) - xubasepanel(i,size(ubase,1)))*cos(twist(i)*pi/180) + ...
        (ubasepanel(i,:) - ubasepanel(i,size(ubase,1)))*sin(twist(i)*pi/180) + xubasepanel(i,size(ubase,1));
    ubasepanel(i,:) = -(xubasepanel(i,:) - xubasepanel(i,size(ubase,1)))*sin(twist(i)*pi/180) + ...
        (ubasepanel(i,:) - ubasepanel(i,size(ubase,1)))*cos(twist(i)*pi/180) + ubasepanel(i,size(ubase,1));
    xlbasepanel(i,:) = (xlbasepanel(i,:) - xlbasepanel(i,size(lbase,1)))*cos(twist(i)*pi/180) + ...
        (lbasepanel(i,:) - lbasepanel(i,size(lbase,1)))*sin(twist(i)*pi/180) + xlbasepanel(i,size(ubase,1));
    lbasepanel(i,:) = -(xlbasepanel(i,:) - xlbasepanel(i,size(lbase,1)))*sin(twist(i)*pi/180) + ...
        (lbasepanel(i,:) - lbasepanel(i,size(lbase,1)))*cos(twist(i)*pi/180) + lbasepanel(i,size(ubase,1));
        %Plotting
%     figure(1);
%     hold on;
%     for j = 1:size(ubase,1)
%        scatter3(xubasepanel(i,j),ubasepanel(i,j),ypanelpos(i),'b.');
%     end
%     for j = 1:size(lbase,1)
%        scatter3(xlbasepanel(i,j),lbasepanel(i,j),ypanelpos(i),'b.');
%     end
%     axis equal tight;
    %Baseline Config Generation end

    %Updating
    
    %Untwisting
    xunewpanel(i,:) = (xubasepanel(i,:) - xubasepanel(i,size(ubase,1)))*cos(-twist(i)*pi/180) + ...
        (ubasepanel(i,:) - ubasepanel(i,size(ubase,1)))*sin(-twist(i)*pi/180) + xubasepanel(i,size(ubase,1));
    unewpanel(i,:) = -(xubasepanel(i,:) - xubasepanel(i,size(ubase,1)))*sin(-twist(i)*pi/180) + ...
        (ubasepanel(i,:) - ubasepanel(i,size(ubase,1)))*cos(-twist(i)*pi/180) + ubasepanel(i,size(ubase,1));
    xlnewpanel(i,:) = (xlbasepanel(i,:) - xlbasepanel(i,size(lbase,1)))*cos(-twist(i)*pi/180) + ...
        (lbasepanel(i,:) - lbasepanel(i,size(lbase,1)))*sin(-twist(i)*pi/180) + xlbasepanel(i,size(ubase,1));
    lnewpanel(i,:) = -(xlbasepanel(i,:) - xlbasepanel(i,size(lbase,1)))*sin(-twist(i)*pi/180) + ...
        (lbasepanel(i,:) - lbasepanel(i,size(lbase,1)))*cos(-twist(i)*pi/180) + lbasepanel(i,size(ubase,1));
    %Translating back
    xunewpanel(i,:) = xunewpanel(i,:) + (tspan-ypanelpos(i))/tan(leangle);
    xlnewpanel(i,:) = xlnewpanel(i,:) + (tspan-ypanelpos(i))/tan(leangle); 
    %Scaling down
    unewpanel(i,:) = unewpanel(i,:)./chordlen(i);
    lnewpanel(i,:) = lnewpanel(i,:)./chordlen(i);
    xunewpanel(i,:) = xunewpanel(i,:)./chordlen(i);
    xlnewpanel(i,:) = xlnewpanel(i,:)./chordlen(i);
    
    % Hicks-Henne Update
    for k = 1:size(ubase)
        sum = 0;
        for j = 1:hh
            m = log(0.5)/log(bump_pos(j));
            sum = sum + dpa(i,j)*(sin(pi*xunewpanel(i,k)^m)^t_b);
            %gradu(i,j) = sin(pi*xunewpanel(i,k)^m)^t_b;
        end
        unewpanel(i,k) = unewpanel(i,k) + sum;
    end
    for k = 1:size(lbase)
        sum = 0;
        for j = 1:hh
            m = log(0.5)/log(bump_pos(j));
            sum = sum + dpb(i,j)*(sin(pi*xlnewpanel(i,k)^m)^t_b);
            %gradu(i,j) = sin(pi*xunewpanel(i,k)^m)^t_b;
        end
        lnewpanel(i,k) = lnewpanel(i,k) + sum;
    end
    
    %Scaling
    unewpanel(i,:) = scale(i).*unewpanel(i,:);
    lnewpanel(i,:) = scale(i).*lnewpanel(i,:);
    xunewpanel(i,:) = scale(i).*xunewpanel(i,:);
    xlnewpanel(i,:) = scale(i).*xlnewpanel(i,:);
    
    %Translating
    xunewpanel(i,:) = xunewpanel(i,:) - dispx(i);
    xlnewpanel(i,:) = xlnewpanel(i,:) - dispx(i);
    unewpanel(i,:) = unewpanel(i,:) + dispy(i);
    lnewpanel(i,:) = lnewpanel(i,:) + dispy(i);
    
    %twisting
    xunewpanel(i,:) = (xunewpanel(i,:) - xunewpanel(i,size(ubase,1)))*cos(twistnew(i)*pi/180) + ...
        (unewpanel(i,:) - unewpanel(i,size(ubase,1)))*sin(twistnew(i)*pi/180) + xunewpanel(i,size(ubase,1));
    unewpanel(i,:) = -(xunewpanel(i,:) - xunewpanel(i,size(ubase,1)))*sin(twistnew(i)*pi/180) + ...
        (unewpanel(i,:) - unewpanel(i,size(ubase,1)))*cos(twistnew(i)*pi/180) + unewpanel(i,size(ubase,1));
    xlnewpanel(i,:) = (xlnewpanel(i,:) - xlnewpanel(i,size(lbase,1)))*cos(twistnew(i)*pi/180) + ...
        (lnewpanel(i,:) - lnewpanel(i,size(lbase,1)))*sin(twistnew(i)*pi/180) + xlnewpanel(i,size(ubase,1));
    lnewpanel(i,:) = -(xlnewpanel(i,:) - xlnewpanel(i,size(lbase,1)))*sin(twistnew(i)*pi/180) + ...
        (lnewpanel(i,:) - lnewpanel(i,size(lbase,1)))*cos(twistnew(i)*pi/180) + lnewpanel(i,size(ubase,1));
    
%     figure(2)
%     hold on;
%     for j = 1:size(ubase,1)
%        scatter3(xunewpanel(i,j),unewpanel(i,j),ypanelpos(i),'r.');
%     end
%     for j = 1:size(lbase,1)
%        scatter3(xlnewpanel(i,j),lnewpanel(i,j),ypanelpos(i),'r.');
%     end
%     axis equal tight;

    
     


end



