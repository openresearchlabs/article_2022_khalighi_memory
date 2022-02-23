%% PanelD
%  Dissimilarity to the initial stable state through time for one randomly 
%  chosen community (community~6) and different memory strengths.
% -------------------------------------------------------------------------
%                   This code solves an N species microbial community model
%                   described by fractional differential equations:
%                   D^mu(Xi)=Xi(bi.Fi-ki.Xi)
%                   where Fi=\prod[Kik^n/(Kik^n+Xk^n)], k=1,...,N and k~=i
%                   D is the fractional Caputo derivative and mu is its order                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs (Except iR, T, T1, and T2, Don't change the input to keep the equilibrium points)
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
%        ------------------------------------------------------------------
%        iR - Number of the random community, e.g. iR=6; The 6th community is choosen out of generated random communities 
%-----------------------------------
% Outputs
%        t - Simulated time interval
%        D - Dissimilarity to the initial stable state through time  
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

t0=0; T=20; % t0 is start time, and T is the end time

load X0(n) % loading equilibrium points based on the parameters

load('Kijfig3e.mat'); % Loading random matrix interaction

Ki=2*ones(N,1); % death rate

n=4; % Hill coefficient

N=15; % number of species

mu=1:-.1:.4; % Order of derivatives

iR=6; % number of the random community

h=0.01; % step size for computing

%% Computing time series of abundances

s = RandStream('twister','Seed',10+iR); % seed for generating random number
    RandStream.setGlobalStream(s)

x0=Xn(:,iR); % equilibrium initial values

RR=randn(N,1);  
b=ones(N,1)+.05*RR; % generating random growth rates

for iM=1:length(mu)

alpha=mu(iM)*ones(N,1); % order of derivatives

[t, x] = FDE_PI12_PC(alpha,@fun3,t0,T,x0,h,[],4);

%%Computing dissimilarity to the initial stable state through time

RelX=x./(ones(N,1)*sum(x)); % relative abundance

X1=RelX(:,1); % initial values
X2=RelX(:,2:end); % Values for comparing with initial values

for it=1:length(t)-1
D(iM,it)=braycd(X1,X2(:,it)); % dissimilarity to the initial values
end

end

%% Plotting

figure 

% Highlight background as perturbation
RX=max(max(D))+0.04; % value for background highlight
f = [1 2 3 4];
v = [T1 0; T2 0; T2 RX; T1 RX];
h1=patch('Faces',f,'Vertices',v,'FaceColor','b', 'EdgeColor','non', 'FaceAlpha',.16);
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');                  

axis tight
xaxiss=1-mu;

startColor = [0, 0, .4];
for ii = 1:length(mu)
    p(ii)=plot(t(2:end), D(ii,:), 'Color', startColor + (ii-.1) / 15);
    hold on;
end

% Settings for plot

ylabel('Dissimilarity to initial state')
xlabel('Time')
set(gca,'FontSize',14)

% Generate legend or title for showing memory
    k =1; ind = 1; 
while k < length(mu)+1
   set(p(ind),'DisplayName',['Memory=',num2str(1-mu(k)),])
   k = k+1; ind = ind+1;
end
hleg1=legend('show');
set(hleg1,'Location','northeast','FontSize',11)
title(['Community ',num2str(iR),])

set(gcf,'renderer','Painters')
