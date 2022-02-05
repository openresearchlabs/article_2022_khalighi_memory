%% Figure S8
clear 
clc

global n N Ki b Kij
%% Inputs
% Coefficients and Conditions

N=2;

order1=1:-.02:.9;
order2=1:-.02:.9;

n=2; % Hill coefficient

Kij=0.1*ones(2); % interaction matrix

Ki=1*ones(N,1); % death rate

T=250; %  final time

b=[1, 2]; % growth rates for cases: False, Pulse, and Periodic

t0=0; % initial time
h=.2; % step size for computing
F=@funGonze;

%%fix points
x2 = 1.9987539067056423719997810984789;
x1=100.*x2.^4 - 200.*x2.^3 + x2.^2 - 2.*x2 + 1;

X0=[x1; x2]; % initial conditions


M1=length(order1);
M2=length(order2);
Resistance=zeros(M1,M2);



for i=1:M1
    for j=1:M2
        tic
        p=34.5; %perturb
while 1
[t,X]=FDE_PI12_PC([order1(i),order2(j)],F,t0,T,X0,h,p);

Df=diff(X(:,end-1:end)');

if Df(1)<=0 && Df(2)>=0
    p=p+0.1;
else
    Resistance(i,j)=p-0.1;
    break
end

end
toc
    end
end

%% plotting
figure

h=heatmap(1-order1,1-order2,Resistance');
h.XLabel = 'Memory of X_B';
h.YLabel = 'Memory of X_R';

ax = gca;
axp = struct(ax);       %you will get a warning
axp.Axes.XAxisLocation = 'top';
