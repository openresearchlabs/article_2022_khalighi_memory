%% Figure 5 (panels b and d) 
%        ------------------------------------------------------------------
%                   This code solves a there species microbial community model
%                   described by fractional differential equations:
%                   D^mu(Xi)=X_i(bi.Fi-ki.Xi)
%                   where Fi=\prod[Kik^n/(Kik^n+Xk^n)], k=1,...,N and k~=i
%                   D is the fractional Caputo derivative and mu is its order  
%
%  For a article titled "Quantifying the impact of ecological memory on the dynamics of interacting communities"         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
global n N b Ki

%% Inputs

load('bb.mat');b=bb; % growth rates 
load('x0.mat'); % initial conditions
load('Kij.mat'); % predifined matrix interactions 

n=2; % Hill coefficient

N=15; % number of species

Blue=1:5;Red=6:10;Green=11:15; % index for separating the species into three groups

Ki=1*ones(N,1); % death rate
t0=0; T=700; % t0 is start time, and T is the end time

Memory=[0;.4;0.148159;0.148160;0.2];
 
h=0.01; % step size for computing

for I=1:length(Memory)

mu=ones(N,1);
mu(Blue)=1-Memory(I); % order of derivatives
% solver for fractional differential equation    
[tt, xx] = FDE_PI12_PC(mu,@fun,t0,T,x0,h);


%% ploting

% defining blindfriendly colors (red and green)
PcG= [0.18,0.40,0.14];
PcR= [0.92,0.27,0.18];

t=tt(1:50:end);
x=xx(:,1:50:end);


%%plotting abundance of species
f=figure;
f.Renderer='painters';

pb=semilogx(t,x(Blue,:),'b');
set(pb,'LineWidth',2.3)
hold on
pr=semilogx(t,x(Red,:),'Color',PcR);
set(pr,'LineWidth',2.3)
pg=semilogx(t,x(Green,:),'Color',PcG);
set(pg,'LineWidth',2.3)
line([t0 T],[0.5,0.5],'LineStyle','--', 'color', 'k')
hold off

% Settings for the general plot
axis([t0 T 0 1])
xlabel('{time}','FontSize',15)

ylabel('Abundnace','FontSize',15)

set(gca,'Fontsize',24)

for i=1:5
    if ismember(i,2:5)==1
set(get(get(pb(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(pr(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(pg(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
% legend({'X_B', 'X_R', 'X_G'},'FontSize',14)


xticks([10^-2 10^0 10^2])
xticklabels({'-2','0','2'})
         
Pos = [795 359 700 457];
set(0, 'DefaultFigurePosition', Pos);
text(.45,.7,['Memory_{B}=',num2str(1-mu(1))],'FontSize',20)
text(.45,.6,['Memory_{R}=',num2str(1-mu(6))],'FontSize',20)
text(.45,.5,['Memory_{G}=',num2str(1-mu(15))],'FontSize',20)

end
