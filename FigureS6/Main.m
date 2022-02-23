%% Figure S6
clear 
clc

global n N Ki b Kij
%% Inputs
%Coefficients and Conditions

N=2;

order1=1:-.02:.84;
order2=1:-.02:.84;

n=2; % Hill coefficient

Kij=0.1*ones(2); % interaction matrix

Ki=1*ones(N,1); % death rate

T=7000; %  final time

b=[1, 2]; % growth rates for cases: False, Pulse, and Periodic

X0=[.9;.2]; % initial conditions

t0=0;
h=.2;
F=@funGonze;

%%fix points
x2 = [0.021687957003156231859766150168713,1.9987539067056423719997810984789];
x1=100.*x2.^4 - 200.*x2.^3 + x2.^2 - 2.*x2 + 1;


M1=length(order1);
M2=length(order2);
ConvergT=zeros(M1,M2);
for i=1:M1
    for j=1:M2
        if j>M2-2
            T=14000;
        end
[t,X]=FDE_PI12_PC([order1(i),order2(j)],F,t0,T,X0,h);

Err=braycd(X(:,end),[x1;x2]);
[~,indFix]=min(Err);

indx=find(braycd(X,[x1(indFix);x2(indFix)])<5e-3);
ConvergT(i,j)=t(indx(1));
    end
end

%% plotting
figure
h=heatmap(1-order1,1-order2,ConvergT');
h.XLabel = 'Memory of X_B';
h.YLabel = 'Memory of X_R';

ax = gca;
axp = struct(ax);       %you will get a warning
axp.Axes.XAxisLocation = 'top';
