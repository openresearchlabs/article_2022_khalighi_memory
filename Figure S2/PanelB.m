%% PanelB
%  Time series of the community 6.
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
%        RelX - Time series of relative abundance   
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

mu=.3; % Order of derivatives
alpha=mu*ones(N,1); 

h=0.01; % step size for computing

iR=6; % number of the random community


%% Computing time series of abundances

s = RandStream('twister','Seed',10+iR); % seed for generating random number
RandStream.setGlobalStream(s)

RR=randn(N,1);
b=ones(N,1)+.05*RR; % generating random growth rates

x0=Xn(:,iR); % equilibrium initial values

[t, x] = FDE_PI12_PC(alpha,@fun3,t0,T,x0,h,[],4);

RelX=x./(ones(N,1)*sum(x)); % relative abundance

%% plotting

% making color plot
Pcolor=zeros(N,3);

p1=hex2rgb('#332288');
   Pcolor(1,:)= p1;
   p2=hex2rgb('#117733');
    Pcolor(2,:)= p2;
    p3=hex2rgb('#44AA99');
        Pcolor(3,:)= p3;
        p4=hex2rgb('#A2C19D');
    Pcolor(4,:)= p4;
    p5=hex2rgb('#DDCC77');
        Pcolor(5,:)= p5;
        p6=hex2rgb('#AF7AF2');
        Pcolor(6,:)= p6;
        p7=hex2rgb('#882255');
        Pcolor(7,:)= p7;
        p8=hex2rgb('#88CCEE');
        Pcolor(8,:)= p8;
        p9=hex2rgb('#F00248');
        Pcolor(9,:)= p9;
        p10=hex2rgb('#A9B222');
        Pcolor(10,:)=p10;
        p11=hex2rgb('#893ED2');
        Pcolor(11,:)= p11;
        p12=hex2rgb('#420F63');
        Pcolor(12,:)= p12;
        p13=hex2rgb('#C0E436');
        Pcolor(13,:)= p13;
        p14=hex2rgb('#DF4298');
    Pcolor(14,:)=  p14;
    p15=hex2rgb('#304B1C');
        Pcolor(15,:)= p15;   

fig=figure;

% set the size of figure
Pos = [300 300 1500 600];
                set(0, 'DefaultFigurePosition', Pos);
                
subplot(1,2,1)
hold on

% Highlight background as perturbation
RX=max(max(RelX))+0.05; % value for background highlight
f = [1 2 3 4];
v = [T1 min(x,[],'all')*.1; T2 min(min(x))*.1; T2 RX; T1 RX];
h1=patch('Faces',f,'Vertices',v,'FaceColor','b', 'EdgeColor','non', 'FaceAlpha',.16);
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');                  
set(gca,'YScale','log')
    
% Plot the relative abundances
p=gobjects(N,1);

for i=1:N
p(i)=semilogy(t(1:20:end),RelX(i,1:20:end)); % relative abundance
p(i).Color=Pcolor(i,:);
end

yticks([10^-4 10^-3 10^-2 10^-1 1])
yticklabels({'-4','-3','-2','-1','0'})

% Settings for plot
set(p,'LineWidth',3)
ylabel('{Log abundance}')
xlabel('{Time}')
set(gca,'FontSize',24)

set(gcf,'renderer','Painters')

%%%%%%%%%%%%%Text the perturbation on plot
text(20,RX,'{3 x b_{J}}','FontSize',23)

text(108,RX-.05,'{J}=','Color','k','FontSize',24)
text(122,RX-.05,'{1}','Color',Pcolor(1,:),'FontSize',24)
text(132,RX-.05,'{2}','Color',Pcolor(2,:),'FontSize',24)
text(142,RX-.05,'{3}','Color',Pcolor(3,:),'FontSize',24)
text(152,RX-.05,'{4}','Color',Pcolor(4,:),'FontSize',24)
text(162,RX-.05,'{5}','Color',Pcolor(5,:),'FontSize',24)
text(172,RX-.05,'{6}','Color',Pcolor(6,:),'FontSize',24)
text(182,RX-.05,'{7}','Color',Pcolor(7,:),'FontSize',24)

hold off

axis tight
xlim([0,T])
% -------------------------------------------------------------------------

% plotting phase plane of final and intial states
subplot(1,2,2)
       set(gca,'FontSize',23)     
hold on
for i=1:N
plot(RelX(i,1),RelX(i,end), '.','Color', Pcolor(i,:), 'MarkerSize',50); % relative abundance
text(RelX(i,1),RelX(i,end),['X_{',num2str(i),'}'],'FontSize',23)
end
% setting for plot
MAX=1.1*max([RelX(:,1);RelX(:,end)]);
line([0,MAX],[0,MAX],'LineStyle','--');
ylabel('Final time', 'FontSize',23)
xlabel('1st steady state', 'FontSize',23)
axis tight
set(gcf,'renderer','Painters')

% Give common xlabel, ylabel and title to your figure
% Create a new axis
ax = axes(fig);
% Specifiy the active side for the axes ax
yyaxis(ax, 'left');
% Specify visibility of the current axis as 'off'
han = gca;
han.Visible = 'off';
% Specify visibility of Title, XLabel, and YLabel as 'on'
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
% Give title, xlabel and left ylabel
% Generate legend or title for showing memory
if all(mu==mu(1))
    tl=title(['Memory=',num2str(1-mu(1))]);
    set(tl,'FontSize',24)
else
    k =1; ind = 1; 
while k < length(mu)+1
   set(p(ind),'DisplayName',['X_{',num2str(k),...
       '}, \mu_{',num2str(k),'}=',num2str(mu(k))])
   k = k+1; ind = ind+1;
end
hleg1=legend('show');
set(hleg1,'Location','northeast','FontSize',11)

indx1=find(t==T1)-1;

X1=x(:,indx1);
X2=x(:,end);

D(iM)=braycd(X1,X2);
toc
end
