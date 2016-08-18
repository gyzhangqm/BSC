%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Hicks-Henne Bump function Parameterization                   %
%                   16 Aug 2016                                      %
%             3D Wing parameterization                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
format long;
load airfoils/oneraup.txt
load airfoils/oneradown.txt
ubase = oneraup(:,2);
xu = oneraup(:,1);
lbase = oneradown(:,2);
xl = oneradown(:,1);

panels = 3;
hh = 5;
t_b = 4;
ypanelpos = zeros(panels,1);
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
twist = zeros(panels,1);
twist(1) = 10;
twist(2) = 5;
twist(3) = 2;


span = 1196.3;
tspan = 2734.5;
xspan = 1578.7;
sweep = 26.7;
leangle = 60;
teangle = 74.219120;
sweep = sweep*pi/180;
leangle = leangle*pi/180;
teangle = teangle*pi/180;
ypanelpos = 0:span/(panels-1):span;
dpa = zeros(panels,hh);
dpb = zeros(panels,hh);
dpa(1,1) = 0.1;
dpa(2,1) = 0.1;
dpa(3,1) = 0.1;
dpb(1,1) = -0.1;
dpb(2,1) = -0.1;
dpb(3,1) = -0.1;

for i = 1:panels
    h = 1/(hh+1);
    bump_pos = 0+h:h:1-h;
    
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
    figure(1);
    hold on;
    for j = 1:size(ubase,1)
       scatter3(xubasepanel(i,j),ubasepanel(i,j),ypanelpos(i),'r*');
    end
    for j = 1:size(lbase,1)
       scatter3(xlbasepanel(i,j),lbasepanel(i,j),ypanelpos(i),'r*');
    end
    axis equal tight;

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
    unewpanel(i,:) = chordlen(i).*unewpanel(i,:);
    lnewpanel(i,:) = chordlen(i).*lnewpanel(i,:);
    xunewpanel(i,:) = chordlen(i).*xunewpanel(i,:);
    xlnewpanel(i,:) = chordlen(i).*xlnewpanel(i,:);
    
    %Translating
    xunewpanel(i,:) = xunewpanel(i,:) - (tspan-ypanelpos(i))/tan(leangle);
    xlnewpanel(i,:) = xlnewpanel(i,:) - (tspan-ypanelpos(i))/tan(leangle);
    
    %twisting
    xunewpanel(i,:) = (xunewpanel(i,:) - xunewpanel(i,size(ubase,1)))*cos(twist(i)*pi/180) + ...
        (unewpanel(i,:) - unewpanel(i,size(ubase,1)))*sin(twist(i)*pi/180) + xunewpanel(i,size(ubase,1));
    unewpanel(i,:) = -(xunewpanel(i,:) - xunewpanel(i,size(ubase,1)))*sin(twist(i)*pi/180) + ...
        (unewpanel(i,:) - unewpanel(i,size(ubase,1)))*cos(twist(i)*pi/180) + unewpanel(i,size(ubase,1));
    xlnewpanel(i,:) = (xlnewpanel(i,:) - xlnewpanel(i,size(lbase,1)))*cos(twist(i)*pi/180) + ...
        (lnewpanel(i,:) - lnewpanel(i,size(lbase,1)))*sin(twist(i)*pi/180) + xlnewpanel(i,size(ubase,1));
    lnewpanel(i,:) = -(xlnewpanel(i,:) - xlnewpanel(i,size(lbase,1)))*sin(twist(i)*pi/180) + ...
        (lnewpanel(i,:) - lnewpanel(i,size(lbase,1)))*cos(twist(i)*pi/180) + lnewpanel(i,size(ubase,1));
    
    figure(2);
    hold on;
    for j = 1:size(ubase,1)
       scatter3(xunewpanel(i,j),unewpanel(i,j),ypanelpos(i),'r*');
    end
    for j = 1:size(lbase,1)
       scatter3(xlnewpanel(i,j),lnewpanel(i,j),ypanelpos(i),'r*');
    end
    axis equal tight;

    
     


end



