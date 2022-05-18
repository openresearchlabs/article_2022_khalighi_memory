%% Figure S11
% CH ER 
clear 
clc
%% Inputs
global A mu

order1=1:-.01:.9; % order of derivatives for CH
order2=1:-.01:.9; % order of derivatives for ER

% initial conditions
ExpRelCH=[0.0077042 0.81112	0.28844	0.58742	0.5202	0.34852	0.51262];
ExpRelER=1-ExpRelCH;
ExpT=0:12:12*length(ExpRelCH)-12;
AbsIniCH=.02*ExpRelCH(1);
AbsIniER=.02*ExpRelER(1);
X0= [AbsIniCH; AbsIniER]; 

mu=[.468 0.151]; % growth rates

t0=0; % initial time
T=510; % final time 
h=.1; % step size for computing
F=@fun; % ODE funcion described by Venturelli et. al. (https://doi.org/10.15252/msb.20178157)
JF=@Jfun; % Jacobian of ODE


A=[-1.242 -.508; 1.191 -1.3219]; % interaction matrix

%% fix points
xx1=[(A(1,2)*mu(2)-A(2,2)*mu(1))/(A(1,1)*A(2,2)-A(1,2)*A(2,1)),...
    (A(2,1)*mu(1)-A(1,1)*mu(2))/(A(1,1)*A(2,2)-A(1,2)*A(2,1))];
xx2=flip(xx1);
xx3=[0,-mu(2)/A(2,2)];
xx4=[-mu(1)/A(1,1),0];

x12=[xx1;xx2;xx3;xx4];

% Jac=JF(1,[x1,x2]);
% eig(Jac);

M1=length(order1);
M2=length(order2);
ConvergT=zeros(M1,M2);
for i=1:M1
    for j=1:M2
[t,X]=FDE_PI2_IM([order1(i),order2(j)],F,JF,t0,T,X0,h);

Err=braycd(X(:,end),x12');
[~,indFix]=min(Err);

indx=find(braycd(X(:,35.7/h:end),x12(indFix,:)')<5e-4);
ConvergT(i,j)=t(indx(1))+35.7;
    end
end

%% Plotting

figure
h=heatmap(1-order1,1-order2,ConvergT');
h.XLabel = 'Memory of CH';
h.YLabel = 'Memory of ER';

ax = gca;
axp = struct(ax);       %you will get a warning
axp.Axes.XAxisLocation = 'top';
