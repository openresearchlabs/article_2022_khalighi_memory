%% PanelC
%  Community dissimilarity (Bray-Curtis) between the start and the end of   
%  the simulation for all ten communities and different memory levels.
% -------------------------------------------------------------------------
%                   This code solves an N species microbial community model
%                   described by fractional differential equations:
%                   D^mu(Xi)=Xi(bi.Fi-ki.Xi)
%                   where Fi=\prod[Kik^n/(Kik^n+Xk^n)], k=1,...,N and k~=i
%                   D is the fractional Caputo derivative and mu is its order                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs (Except T, T1, and T2, Don't change the input to keep the equilibrium points)
%
%        mu - Order of derivatives, e.g. mu=0.7*ones(1,N);  % 1-Memory
%        ------------------------------------------------------------------
%        n -  Hill coefficient, e.g. n=4;                              
%        ------------------------------------------------------------------
%        N -  Number of Species, N=15;
%        ------------------------------------------------------------------
%        Ki - Death rate, e.g. Ki=1*ones(N,1); 
%        ------------------------------------------------------------------
%        T - Final time, e.g. T=200;
%        ------------------------------------------------------------------
%        T1 - When the pulse perturbation starts 
%        T2 - When the pulse perturbation ends
%
%-----------------------------------
% Outputs
%        t - Simulated time interval
%        D - Bray-Curtis between the start and the end of the simulation   
%-----------------------------------
% Moein Khalighi - April 2021
%
%
%  Please, report any problem or comment to :
%          moein dot khalighi at utu dot fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
global n N Ki mu T1 T2 b

%% Inputs

T1=20; % When the pulse perturbation starts 
T2=70; % When the pulse perturbation ends

t0=0; T=200; % t0 is start time, and T is the end time

load X0(n) % loading equilibrium points based on the parameters

load('Kijfig3e.mat'); % Loading random matrix interaction

Ki=2*ones(N,1); % death rate

n=4; % Hill coefficient

N=15; % number of species

h=0.01; % step size for computing

mu=1:-.1:.4; % Order of derivatives

%%
for iR=1:10
    s = RandStream('twister','Seed',10+iR); % seed for generating random number
    RandStream.setGlobalStream(s)

RR=randn(N,1);  
b=ones(N,1)+.05*RR; % generating random growth rates

x0=Xn(:,iR); % equilibrium initial values

for iM=1:length(mu)

alpha=mu(iM)*ones(N,1); % order of derivatives

[t, x] = FDE_PI12_PC(alpha,@fun3,t0,T,x0,h,[],4);%%Computing abundance for each random community

RelX=x./(ones(N,1)*sum(x));% relative abundance

% Computing Bray-Curtis between the start and the end of the simulation 
indx1=find(t==T1)-1;

X1=RelX(:,indx1);
X2=RelX(:,end);

D(iR,iM)=braycd(X1,X2);

end
end
%% Plotting

xaxiss=1-mu; 

boxplot(D,xaxiss)

% Settings for plot

ylabel('{Dissimilarity between \newline initial and final states}')
xlabel('{Memory}')
set(gca,'FontSize',14)
