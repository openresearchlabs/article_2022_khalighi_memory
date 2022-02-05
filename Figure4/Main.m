%% Figure 4
%        ------------------------------------------------------------------
%                   This code solves a there species microbial community model
%                   described by fractional differential equations:
%                   D^mu(Xi)=X_i(bi.Fi-ki.Xi)
%                   where Fi=\prod[Kik^n/(Kik^n+Xk^n)], k=1,...,N and k~=i
%                   D is the fractional Caputo derivative and mu is its order  
%
%  For a article titled "Quantifying the impact of ecological memory on the dynamics of interacting communities"         
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
%        b - Growth rates (including stochastic perturbations), e.g. b=[1, .95, 1.05];
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
global n N Ki b Kij T x0 bb
%% Coefficients and Conditions

%memory values
Memory=[0;.1;.2;.3;.4;0.0842041;0.0842046;0.0842048]; 
mu=ones(1,3)-Memory; % order of derivatives 

n=2; % Hill coefficient

N=3; % number of species

Kij=0.1*ones(N); % interaction matrix

Ki=1*ones(N,1); % death rate

T=700; %  final time

x0=1/3*[1,1,1]'; % initial conditions

Fun=@fun3;
        clear b        
        load('b3OUP.mat'); % growth rates are generated under an oup process
        bb=@(t,N)b(t,N);
        B=b;
        
        t0=0; % initial time
h=0.01; % step size for computing

% solver for fractional differential equation
for i=1:length(Memory)

[t, x] = FDE_PI12_PC(mu(i,:),Fun,t0,T,x0,h);
X(i,:,:)=x;

end

%% plotting
x=zeros(3,length(t));
for i=1:length(Memory)
x(:,:)=X(i,:,:);
RelX=x./(ones(N,1)*sum(x)); % Relative abundances

f=figure;
f.Renderer='painters';

% defining blindfriendly colors (red and green)
PcR= [0.92,0.27,0.18];
PcG= [0.18,0.40,0.14];
          
    %%plotting relative abundance of species
       
            p=plot(t,RelX(1,:),'b',t,RelX(2,:),'r',t,RelX(3,:),'g');
       ylabel('Relative abundance')


% Settings for the general plot
set(p,'LineWidth',5)
set(gca,'FontSize',23, 'FontWeight', 'bold')
xlabel('Time')

% % Generate legend for showing memory of each species
legend(['X_{B}/X_{total}, Memory_B=',num2str(1-mu(i,1))],...
    ['X_{R}/X_{total}, Memory_R=',num2str(1-mu(i,2))],...
    ['X_{G}/X_{total}, Memory_G=',num2str(1-mu(i,3))],'Location', 'Best');

axis([0,700,0,1])
end

  %%plot of growth rate with stochastic perturbation

f=figure;
f.Renderer='painters';
    
    tt=1:T; % time simulating with 1 step size  
    
    % Plot the simulated growth rates
    
    pb=plot(tt,B(1:T,1),'b',tt,B(1:T,2),'r',tt,B(1:T,3),'g');
    
    % Settings for plot
    set(pb,'LineWidth',4)
    ylabel('Growth Rates: OUP')
    xlabel('Time')
    legend('b_B','b_R', 'b_G')
    set(gca,'FontSize',23, 'FontWeight', 'bold')
  