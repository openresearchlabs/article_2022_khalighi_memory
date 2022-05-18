%% Figure S5 panel e
% CH ER 
clear 
clc

global A mu

%% Inputs

p=1:-.1:.6; % order of derivatives

ExpRelCH=[0.0077042 0.81112	0.28844	0.58742	0.5202	0.34852	0.51262];
% ExpRelCH=[0.67341	0.94866	0.38006	0.65739	0.53902	0.42119	0.63194];
% ExpRelCH=[0.32481 0.46857 0.42027 0.55545 0.42054 0.53746 0.42026];
ExpRelER=1-ExpRelCH;
ExpT=0:12:12*length(ExpRelCH)-12;

% initial conditions
AbsIniCH=.02*ExpRelCH(1);
AbsIniER=.02*ExpRelER(1);
X0= [AbsIniCH; AbsIniER]; 
mu=[.468 0.151]; % growth rates

t0=0; % intial time
T=100; % final time
h=.1; % step size for computing
F=@fun; % ODE funcion described by Venturelli et. al. (https://doi.org/10.15252/msb.20178157)
JF=@Jfun; % Jacobian of ODE


A=[-1.242 -.508; 1.191 -1.3219]; % interaction coefficients


%% Solver for fractional differential equation 
for j=1:length(p)
    order=p(j)*ones(2,1);
[t,x]=FDE_PI2_IM(order,F,JF,t0,T,X0,h);
X(j,:,:)=x;

end

%% plotting
f1=figure;
f1.Renderer='painters';

% relative abundances
x1(:,:)=X(1,:,:)./sum(X(1,:,:));
x2(:,:)=X(2,:,:)./sum(X(2,:,:));
x3(:,:)=X(3,:,:)./sum(X(3,:,:));
x4(:,:)=X(4,:,:)./sum(X(4,:,:));
x5(:,:)=X(5,:,:)./sum(X(5,:,:));

hold on
p1=plot(t,x1(1,:),'Color',[0,0,.7],'LineWidth',3,'DisplayName','0');
p2=plot(t,x2(1,:),'Color',[0,0,.9],'LineWidth',3,'DisplayName','0.1');
p3=plot(t,x3(1,:),'Color',[0,0.5,1],'LineWidth',3,'DisplayName','0.2');
p4=plot(t,x4(1,:),'Color',[0,.7,1],'LineWidth',3,'DisplayName','0.3');
p5=plot(t,x5(1,:),'Color',[0,.9,1],'LineWidth',3,'DisplayName','0.4');
pEx=plot(ExpT,ExpRelCH,'sb','LineWidth',2,'MarkerSize',10,'DisplayName','Exp');

p6=plot(t,x1(2,:),'Color',[.4,0,0],'LineWidth',3,'DisplayName','0');
p7=plot(t,x2(2,:),'Color',[.6,0,0],'LineWidth',3,'DisplayName','0.1');
p8=plot(t,x3(2,:),'Color',[.8,0,0],'LineWidth',3,'DisplayName','0.2');
p9=plot(t,x4(2,:),'Color',[1,0,0],'LineWidth',3,'DisplayName','0.3'); 
p10=plot(t,x5(2,:),'Color',[1,.3,0],'LineWidth',3,'DisplayName','0.4');
pEx1=plot(ExpT,ExpRelER,'sr','LineWidth',2,'MarkerSize',10,'DisplayName','Exp');


set(gca,'FontSize',14)



xlabel('Time')
ylabel('Abundance')

set(gca,'FontSize',14)
% legend('0','0.1','0.2', '0.3', '0.4')

leg = legend('show');
title(leg,'CH    \mu   ER   \mu ')
leg.NumColumns = 2;
leg.Location= 'best';
