%% Figure 6b 
% This code solve a logistic growth curve with memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%% Inputs
F=@fun1; % logistic growth model 
mu=1:-.1:.6; % order of derivative
h=0.01; % step size for comuting
t0=0; % initial time
T=10; %  final time
x0=.1; % initial condition
 
%% solver for fractional differential equation
for i=1:length(mu)
    [t,x] = FDE_PI12_PC(mu(i),F,t0,T,x0,h);
    X(i,:)=x;
end

%% plotting
figure
p1=plot(t,X);
xlabel('Time')
ylabel('Logistic growth curve')
set(p1,'LineWidth',3)
set(gca,'FontSize',14)
legend('0','0.1','0.2', '0.3', '0.4')

leg = legend('show');
title(leg,'Memory')
