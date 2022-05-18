%% Figure S5 panel a and b
%        ------------------------------------------------------------------
%                   This code solves a there species microbial community model
%                   described by fractional differential equations:
%                   D^mu(Xi)=X_i(bi.Fi-ki.Xi)
%                   where Fi=\prod[Kik^n/(Kik^n+Xk^n)], k=1,...,N and k~=i
%                   D is the fractional Caputo derivative and mu is its order  
%
%  For a article titled "Quantifying the impact of ecological memory on the dynamics of interacting communities"         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs                   
%        ------------------------------------------------------------------
%        mu - Order of derivatives, [mu_B,mu_R,mu_G]  0<mu(i)=<1, e.g. mu=[1,.2,1];
%        ------------------------------------------------------------------
%        n -  Hill coefficient, e.g. n=2;
%        ------------------------------------------------------------------
%        N -  Number of Species, e.g. N=3;
%        ------------------------------------------------------------------
%        Kij - Interation matrix, e.g. Kij=0.1*ones(N);
%        ------------------------------------------------------------------
%        Ki - Death rate, e.g. Ki=1*ones(N,1);
%        ------------------------------------------------------------------
%        T - Final time, e.g. T=600;
%        ------------------------------------------------------------------
%        x0 - Initial conditions, e.g. x0=[1/3;1/3;1/3];
%        ------------------------------------------------------------------
%        b - Growth rates (including pulse perturbations), e.g. b=[1, .95, 1.05];
%
%---------------------------------------
% Outputs
%        t - Simulated time interval
%        x - Species abundances 

clear 
clc

global n N Ki b Kij
%% Inputs
% Coefficients and Conditions

N=2;

p=1:-.1:.6; % order of derivatives

n=2; % Hill coefficient

Kij=0.1*ones(2); % interaction matrix

Ki=1*ones(N,1); % death rate

T=100; %  final time

b=[1, 2]; % growth rates

% initial conditions
X0=[0.9;0.15]; % for panel a 
% X0=[.8;.2]; % for panel b

t0=0; % initial time
h=0.2; % step size for computing
F=@funGonze; %ODE model

%% solver for fractional differential equation

for j=1:length(p)
    order=p(j)*ones(2,1);
[t,x]=FDE_PI12_PC(order,F,t0,T,X0,h);
X(j,:,:)=x;
end

%% plotting

f1=figure;
f1.Renderer='painters';

% relative abundances
x1(:,:)=X(1,:,:)./sum(X(1,:,:));x1=x1(:,1:10:end);
x2(:,:)=X(2,:,:)./sum(X(2,:,:));x2=x2(:,1:10:end);
x3(:,:)=X(3,:,:)./sum(X(3,:,:));x3=x3(:,1:10:end);
x4(:,:)=X(4,:,:)./sum(X(4,:,:));x4=x4(:,1:10:end);
x5(:,:)=X(5,:,:)./sum(X(5,:,:));x5=x5(:,1:10:end);
t=t(:,1:10:end);

hold on
p1=plot(t,x1(1,:),'Color',[0,0,.7],'LineWidth',3,'DisplayName','0');
p2=plot(t,x2(1,:),'Color',[0,0,.9],'LineWidth',3,'DisplayName','0.1');
p3=plot(t,x3(1,:),'Color',[0,0.5,1],'LineWidth',3,'DisplayName','0.2');
p4=plot(t,x4(1,:),'Color',[0,.7,1],'LineWidth',3,'DisplayName','0.3');
p5=plot(t,x5(1,:),'Color',[0,.9,1],'LineWidth',3,'DisplayName','0.4');

p6=plot(t,x1(2,:),'Color',[.4,0,0],'LineWidth',3,'DisplayName','0');
p7=plot(t,x2(2,:),'Color',[.6,0,0],'LineWidth',3,'DisplayName','0.1');
p8=plot(t,x3(2,:),'Color',[.8,0,0],'LineWidth',3,'DisplayName','0.2');
p9=plot(t,x4(2,:),'Color',[1,0,0],'LineWidth',3,'DisplayName','0.3'); 
p10=plot(t,x5(2,:),'Color',[1,.3,0],'LineWidth',3,'DisplayName','0.4');

set(gca,'FontSize',14)



xlabel('Time')
ylabel('Abundance')

set(gca,'FontSize',14)
% legend('0','0.1','0.2', '0.3', '0.4')

leg = legend('show');
title(leg,'X1    \mu   X2   \mu ')
leg.NumColumns = 2;
