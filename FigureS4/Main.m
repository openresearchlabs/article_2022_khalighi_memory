%% Figure S4
%        ------------------------------------------------------------------
%                   This code solves a there species microbial community model
%                   described by fractional differential equations:
%                   D^mu(Xi)=X_i(bi.Fi-ki.Xi)
%                   where Fi=\prod[Kik^n/(Kik^n+Xk^n)], k=1,...,N and k~=i
%                   D is the fractional Caputo derivative and mu is its order  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs                   
%        ------------------------------------------------------------------
%        mu - Order of derivatives, [mu_B,mu_R,mu_G]  0<mu(i)=<1, e.g. mu=[1,.2,1];
%        ------------------------------------------------------------------
%        n -  Hill coefficient, e.g. n=2;
%        ------------------------------------------------------------------
%        N -  Number of Species, e.g. N=3;
%        ------------------------------------------------------------------
%        Kij - Interation matrix, e.g. Kij=0.1*ones(N);
%        ------------------------------------------------------------------
%        Ki - Death rate, e.g. Ki=1*ones(N,1);
%        ------------------------------------------------------------------
%        T - Final time, e.g. T=600;
%        ------------------------------------------------------------------
%        x0 - Initial conditions, e.g. x0=[1/3;1/3;1/3];
%        ------------------------------------------------------------------
%        Perturbation - This is changes of the species growth rates, e.g. Perturbation='OUP';
%                       Possible usages 'False', 'Pulse, 'Periodic', 'OUP', and 'OUP_new'                      
%                       False: No perturbation
%                       Pulse1: A pulse in (20,60) for Figure 2a
%                       Pulse2: Similar to Pulse1 with a greater strength for Figure 2b
%                       Pulse3: Two pulses in (60,100) and (200,330) for Figure 3a
%                       Pulse4: Two pulses in (60,100) and (400,530) for Figure S2a
%                       Periodic: Periodic perturbation with 20 span for Figure 3b
%                       OUP: Stochastic pertubation used in the paper; requirement: T=<700, For Figure 4b-c & S2b                          
%                       OUP_new: New generating stochastic perturbation
%        ------------------------------------------------------------------
%        b - Growth rates for cases: False, Pulse, and Periodic, e.g. b=[1, .95, 1.05];
%
%---------------------------------------
% Outputs
%        t - Simulated time interval
%        x - Species abundances 
%        B - Growth rates including perturbation
%
%
%  Please, report any problem or comment to :
%          moein dot khalighi at utu dot fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
clc
global n N Ki b Kij 
%% Inputs
% Coefficients and Conditions

% [mu_B,mu_R,mu_G] Order of derivatives,  0<mu(i)=<1
mu1=[1,1,1-.09105];  % for stochastic (panel b)
mu2=[1,1,1-.2];      % for two pulses (panel a)

n=2; % Hill coefficient

N=3; % number of species

Kij=0.1*ones(N); % interaction matrix

Ki=1*ones(N,1); % death rate

T=700; %  final time

Perturbation='OUP'; % Possible usages 'False', 'Pulse1', 'Pulse2', 'Pulse3', 'Pulse4', 'Periodic', 'OUP', and 'OUP_new'

b=[1, .95, 1.05]; % growth rates for cases: False, Pulse, and Periodic

x01=[.99,.01,.01]'; % initial conditions
x02=1/3*[1,1,1]'; % initial conditions

Fun1=@fun14;
Fun2=@fun3;
        clear b        
        load('b3OUP.mat');
        bb=@(t,N)b(t,N);
        B=b;
        
        t0=0; % initial time
h=0.01; % step size for computing


%% solver for fractional differential equation
[t, x1] = FDE_PI12_PC(mu1,Fun1,t0,T,x01,h);
[~, x2] = FDE_PI12_PC(mu2,Fun2,t0,T,x02,h);

%% plotting

   %%plotting relative abundance of species
    figure
     f = [1 2 3 4];
    v = [60 0; 100 0; 100 1; 60 1];
    v1 = [400 0; 530 0; 530 1; 400 1];
    h1=patch('Faces',f,'Vertices',v,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
    h2=patch('Faces',f,'Vertices',v1,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.36);
    
    RelX1=x1./(ones(N,1)*sum(x1)); % Relative abundances
    
    % Plot the relative abundances
    p=plot(t,RelX1(1,:),'b',t,RelX1(2,:),'r',t,RelX1(3,:),'g');
    
    % Remove highlighted bars from the legend
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
       ylabel('Relative abundance')
xlabel('Time')

set(p,'LineWidth',5)
set(gca,'FontSize',23, 'FontWeight', 'bold')
xlabel('Time')
%        ------------------------------------------------------------------

%%
RelX2=x2./(ones(N,1)*sum(x2)); % Relative abundances

f=figure;
f.Renderer='painters';

% defining blindfriendly colors (red and green)
PcR= [0.92,0.27,0.18];
PcG= [0.18,0.40,0.14];
     
   %%plotting relative abundance of species
       
            p=plot(t,RelX2(1,:),'b',t,RelX2(2,:),'r',t,RelX2(3,:),'g');
       ylabel('Relative abundance')


% Settings for the general plot
set(p,'LineWidth',5)
set(gca,'FontSize',23, 'FontWeight', 'bold')
xlabel('Time')

% % Generate legend for showing memory of each species
legend(['X_{B}/X_{total}, Memory_B=',num2str(1-mu2(1))],...
    ['X_{R}/X_{total}, Memory_R=',num2str(1-mu2(2))],...
    ['X_{G}/X_{total}, Memory_G=',num2str(1-mu2(3))],'Location', 'Best');





