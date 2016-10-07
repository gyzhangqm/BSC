clear all;
close all;
clc;
format longe;

load sequential.txt
load plotdata.txt

x = sequential(:,1);
dp1 = sequential(:,2);
dp2 = sequential(:,3);
dp3 = sequential(:,4);
dp4 = sequential(:,5);
dp5 = sequential(:,6);
dp6 = sequential(:,7);
dp7 = sequential(:,8);
dp8 = sequential(:,9);
dp9 = sequential(:,10);
dp10 = sequential(:,11);
dp11 = sequential(:,12);
dp12 = sequential(:,13);
dp13 = sequential(:,14);
obj = sequential(:,15);

figure (1);
hold on;
plot(x,dp1,'o-');
plot(x,dp2,'o-');
plot(x,dp3,'o-');
plot(x,dp4,'o-');
plot(x,dp5,'o-');
plot(x,dp6,'o-');
plot(x,dp7,'o-');
plot(x,dp8,'o-');
plot(x,dp9,'o-');
plot(x,dp10,'o-');
plot(x,dp11,'o-');
plot(x,dp12,'o-');
legend('var1','var2','var3','var4','var5','var6','var7','var8','var9','var10','var11','var12');

axis([1,25,-0.02,0.0125]);
xlabel('Iterations \rightarrow','FontSize',24,'FontWeight','bold','Color','k');
ylabel('\rightarrow','FontSize',24,'FontWeight','bold','Color','k');
title('Variation of Design variables with each iteration','FontSize',24,'FontWeight','bold','Color','k');

figure (2);
plot(x,obj,'o-');
xlabel('Iterations \rightarrow','FontSize',24,'FontWeight','bold','Color','k');
ylabel('\rightarrow','FontSize',24,'FontWeight','bold','Color','k');
title('Variation of Objective function with each iteration','FontSize',24,'FontWeight','bold','Color','k');

