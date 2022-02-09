%% PanelA
%  Dissimilarity to the initial stable state through time for all ten
%  communities, for three different memory strengths.
%  We simulated ten communities of 15 species, each withrandom  interaction  matrices.
%  A similar level of commensurate memory is applied to all ten communities. 
%  Every community is initially in a stable state of the system, and a perturbation 
%  is imposed by multiplying the growth rates of half of the species (b1, ... , b7) by 3.
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

t0=0; T=200; % t0 is start time, and T is the end time

load X0(n) % loading equilibrium points based on the parameters

load('Kijfig3e.mat'); % Loading random matrix interaction

Ki=2*ones(N,1); % death rate

n=4; % Hill coefficient

N=15; % number of species

mu=.7; % Order of derivatives
alpha=mu*ones(N,1); 

h=0.01; % step size for computing

%% Computing

for iR=1:10
    s = RandStream('twister','Seed',10+iR); % seed for generating random number
    RandStream.setGlobalStream(s)

RR=randn(N,1);
b=ones(N,1)+.05*RR; % generating random growth rates

x0=Xn(:,iR); % equilibrium initial values

%%Computing abundance for each random community

[t, x] = FDE_PI12_PC(alpha,@fun3,t0,T,x0,h,[],4);

%%Computing dissimilarity to the initial stable state through time

RelX=x./(ones(N,1)*sum(x)); % relative abundance

X1=RelX(:,1); % initial values
X2=RelX(:,2:end); % Values for comparing with initial values

for it=1:length(t)-1
D(iR,it)=braycd(X1,X2(:,it)); % dissimilarity to the initial values
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
hold on

p=plot(t(2:end),D);
set(p,'LineWidth',2)

% Settings for plot

ylabel('{Dissimilarity to initial state}')
xlabel('{Time}')
set(gca,'FontSize',14)

% Generate legend or title for showing memory
    k =1; ind = 1; 
while k < iR+1
   set(p(ind),'DisplayName',['Community',num2str(k),])
   k = k+1; ind = ind+1;
end
hleg1=legend('show');
set(hleg1,'Location','northeast','FontSize',11)
title(['Memory=',num2str(1-mu),])

set(gcf,'renderer','Painters')
