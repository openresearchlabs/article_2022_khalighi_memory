%% Figure S5 panel f
% BT CH
clear 
clc

global A mu

%% Inputs

p=1:-.1:.6; %order of derivatives

ExpRelBT=[0.17592	0.70232	0.74131	0.9864	0.98978	0.99893	0.99898];
ExpRelCH=1-ExpRelBT;
ExpT=0:12:12*length(ExpRelBT)-12;

% initial conditions
AbsIniBT=.02*ExpRelBT(1);
AbsIniCH=.02*ExpRelCH(1);
X0= [AbsIniBT; AbsIniCH];
mu=[.626 0.468]; % growth rates

t0=0; % initial time
T=100; % final time 
h=.1; % step size for computing
F=@fun; % ODE funcion described by Venturelli et. al. (https://doi.org/10.15252/msb.20178157)
JF=@Jfun; % Jacobian of ODE

A=[-0.9597 -0.0727; -0.5906 -1.242]; % interaction coefficients


%% Solver for fractional differential equation

for j=1:length(p)
    order=p(j)*ones(2,1);
[t,x]=FDE_PI2_IM(order,F,JF,t0,T,X0,h);
X(j,:,:)=x;

end

%% plotting

f1=figure;
f1.Renderer='painters';

x1(:,:)=X(1,:,:)./sum(X(1,:,:));x1=x1(:,1:20:end);
x2(:,:)=X(2,:,:)./sum(X(2,:,:));x2=x2(:,1:20:end);
x3(:,:)=X(3,:,:)./sum(X(3,:,:));x3=x3(:,1:20:end);
x4(:,:)=X(4,:,:)./sum(X(4,:,:));x4=x4(:,1:20:end);
x5(:,:)=X(5,:,:)./sum(X(5,:,:));x5=x5(:,1:20:end);
t=t(:,1:20:end);

hold on
p1=plot(t,x1(1,:),'Color',[0,0,.7],'LineWidth',3,'DisplayName','0');
p2=plot(t,x2(1,:),'Color',[0,0,.9],'LineWidth',3,'DisplayName','0.1');
p3=plot(t,x3(1,:),'Color',[0,0.5,1],'LineWidth',3,'DisplayName','0.2');
p4=plot(t,x4(1,:),'Color',[0,.7,1],'LineWidth',3,'DisplayName','0.3');
p5=plot(t,x5(1,:),'Color',[0,.9,1],'LineWidth',3,'DisplayName','0.4');
pEx=plot(ExpT,ExpRelBT,'sb','LineWidth',2,'MarkerSize',10,'DisplayName','Exp');

p6=plot(t,x1(2,:),'Color',[.4,0,0],'LineWidth',3,'DisplayName','0');
p7=plot(t,x2(2,:),'Color',[.6,0,0],'LineWidth',3,'DisplayName','0.1');
p8=plot(t,x3(2,:),'Color',[.8,0,0],'LineWidth',3,'DisplayName','0.2');
p9=plot(t,x4(2,:),'Color',[1,0,0],'LineWidth',3,'DisplayName','0.3'); 
p10=plot(t,x5(2,:),'Color',[1,.3,0],'LineWidth',3,'DisplayName','0.4');
pEx2=plot(ExpT,ExpRelCH,'sr','LineWidth',2,'MarkerSize',10,'DisplayName','Exp');

set(gca,'FontSize',14)



xlabel('Time')
ylabel('Abundance')

set(gca,'FontSize',14)
% legend('0','0.1','0.2', '0.3', '0.4')

leg = legend('show');
title(leg,'BT    \mu   CH   \mu ')
leg.NumColumns = 2;
