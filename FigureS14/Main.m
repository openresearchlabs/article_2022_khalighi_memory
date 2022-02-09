%% Figure S14
clear 
clc

global A mu
%% inputs
N=2;

order1=1:-.02:.9;% order of derivatives for BU
order2=1:-.02:.9; % order of derivatives for BT

mu=[0.599 0.626];%growth rates

t0=0; % initial time

h=.1;% step size for computing
F=@fun; % ODE funcion described by Venturelli et. al. (https://doi.org/10.15252/msb.20178157)
JF=@Jfun; % Jacobian of ODE

A=[-0.9059 -0.9377;-0.972 -0.9597]; % interaction coefficients

T=1500; %  final time

%% fix points
xx1=[(A(1,2)*mu(2)-A(2,2)*mu(1))/(A(1,1)*A(2,2)-A(1,2)*A(2,1));...
    (A(2,1)*mu(1)-A(1,1)*mu(2))/(A(1,1)*A(2,2)-A(1,2)*A(2,1))];
xx3=[1e-3;-mu(2)/A(2,2)];
X0=xx3; % initial conditions
p=.2;

M1=length(order1);
M2=length(order2);
ConvergT=zeros(M1,M2);

IndxP=80/h;

for i=1:M1
    for j=1:M2
        
[t,X]=FDE_PI2_IM([order1(i),order2(j)],F,JF,t0,T,X0,h,p);

% indx=find(braycd(X(:,IndxP:end),X0)<7.7e-4);
indx=find(braycd(X(:,IndxP:end),X0)<2e-2);
ConvergT(i,j)=t(indx(1))+80;
    end
end

%% plotting
figure
h=heatmap(1-order1,1-order2,ConvergT');
h.XLabel = 'Memory of BU';
h.YLabel = 'Memory of BT';

ax = gca;
axp = struct(ax);       %you will get a warning
axp.Axes.XAxisLocation = 'top';
