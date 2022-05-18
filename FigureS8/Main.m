%% Figure S8
clear 
clc

global n N Ki b Kij
%% Inputs
% Coefficients and Conditions


%convergence interval
Tol= 2e-2; % for panel a
% Tol= 1e-6; % for panel b


N=2;

order1=1:-.02:.9;
order2=1:-.02:.9;

n=2; % Hill coefficient

Kij=0.1*ones(2); % interaction matrix

Ki=1*ones(N,1); % death rate

T=500; %  final time

b=[1, 2]; % growth rates for cases: False, Pulse, and Periodic

t0=0;
h=.1;
F=@funGonze; %ODE function model 1

%%fix points
x2 = 1.9987539067056423719997810984789;
x1=100.*x2.^4 - 200.*x2.^3 + x2.^2 - 2.*x2 + 1;

X0=[x1; x2]; % initial conditions
p=34;

M1=length(order1);
M2=length(order2);
ConvergT=zeros(M1,M2);

IndxP=100/h;

for i=1:M1
    for j=1:M2
        
[t,X]=FDE_PI12_PC([order1(i),order2(j)],F,t0,T,X0,h,p);

indx=find(braycd(X(:,IndxP:end),X0)<Tol);
ConvergT(i,j)=t(indx(1))+100;
    end
end

%% Plotting
figure

h=heatmap(1-order1,1-order2,ConvergT');
h.XLabel = 'Memory of X_B';
h.YLabel = 'Memory of X_R';

ax = gca;
axp = struct(ax);       %you will get a warning
axp.Axes.XAxisLocation = 'top';
