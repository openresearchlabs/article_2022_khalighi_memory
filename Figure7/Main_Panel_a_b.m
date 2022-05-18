%% Figure 7 panel A and B
% BU BT  
clear 
clc

global A mu

%% Inputs

p=1:-.05:.8; %order of derivatives

ExpRelBT=[0.5873 0.62714 0.63223 0.66755 0.67419 0.70986];
ExpRelBU=1-ExpRelBT;
ExpT=0:12:12*length(ExpRelBT)-12;

% initial conditions
AbsIniBT=.02*ExpRelBT(1);
AbsIniBU=.02*ExpRelBU(1);
% X0= [AbsIniBU; AbsIniBT]; % for panel A
X0= [AbsIniBU+.005; AbsIniBT]; % for panel B
mu=[0.599 0.626]; % growth rates

t0=0; % Intial time
T=1000; % Final time
h=.1; % step size for computing
F=@fun; % ODE funcion described by Venturelli et. al. (https://doi.org/10.15252/msb.20178157)
JF=@Jfun; % Jacobian of ODE

A=[-0.9059 -0.9377;-0.972 -0.9597]; % interaction coefficients

%% solver for fractional differential equation

for j=1:length(p)
    order=p(j)*ones(2,1);
[t,x]=FDE_PI2_IM(order,F,JF,t0,T,X0,h);
X(j,:,:)=x;
end


%% plotting
f1=figure;
f1.Renderer='painters';

%relative abundances
x1(:,:)=X(1,:,:)./sum(X(1,:,:));x1=x1(:,1:20:end);
x2(:,:)=X(2,:,:)./sum(X(2,:,:));x2=x2(:,1:20:end);
x3(:,:)=X(3,:,:)./sum(X(3,:,:));x3=x3(:,1:20:end);
x4(:,:)=X(4,:,:)./sum(X(4,:,:));x4=x4(:,1:20:end);
x5(:,:)=X(5,:,:)./sum(X(5,:,:));x5=x5(:,1:20:end);
t=t(:,1:20:end);
hold on
p1=plot(t,x1(1,:),'Color',[0,0,.7],'LineWidth',2,'DisplayName','0');
p2=plot(t,x2(1,:),'Color',[0,0,.9],'LineWidth',2,'DisplayName','0.05');
p3=plot(t,x3(1,:),'Color',[0,0.5,1],'LineWidth',2,'DisplayName','0.1');
p4=plot(t,x4(1,:),'Color',[0,.7,1],'LineWidth',2,'DisplayName','0.15');
p5=plot(t,x5(1,:),'Color',[0,.9,1],'LineWidth',2,'DisplayName','0.2');
% pEx=plot(ExpT,ExpRelBU,'sb','LineWidth',2,'MarkerSize',10,'DisplayName','Exp');

p6=plot(t,x1(2,:),'Color',[.4,0,0],'LineWidth',2,'DisplayName','0');
p7=plot(t,x2(2,:),'Color',[.6,0,0],'LineWidth',2,'DisplayName','0.05');
p8=plot(t,x3(2,:),'Color',[.8,0,0],'LineWidth',2,'DisplayName','0.1');
p9=plot(t,x4(2,:),'Color',[1,0,0],'LineWidth',2,'DisplayName','0.15'); 
p10=plot(t,x5(2,:),'Color',[1,.3,0],'LineWidth',2,'DisplayName','0.2');
% pEx2=plot(ExpT,ExpRelBT,'sr','LineWidth',2,'MarkerSize',10,'DisplayName','Exp');

set(gca,'FontSize',14)


% axis tight
% axis([550,600, 0,1])
axis([950,1000, 0,1])

% axis([0,T, 0,1])
xlabel('Time (hr)')
% ylabel('Abundance')
set(gca,'ytick',[])
yyaxis right

set(gca,'FontSize',14)
% legend('0','0.1','0.2', '0.3', '0.4')

leg = legend('show');
title(leg,'BU     \mu    BT    \mu   ')
leg.NumColumns = 2;
