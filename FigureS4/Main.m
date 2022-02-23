%% Figure S4
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
%        b - Growth rates (including pulse perturbations), e.g. b=[1, .95, 1.05];
%
%---------------------------------------
% Outputs
%        t - Simulated time interval
%        x - Species abundances 
%        B - Growth rates including perturbation
%---------------------------------------
%
%
%  Please, report any problem or comment to :
%          moein dot khalighi at utu dot fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
clc
global n N b Ki Kij

%% Inputs

XB0=.005:.01:.05; % Initial abundances for blue species 
NXB0=length(XB0);
n=2; % Hill coefficient
N=3; % number of species
Kij=0.1*ones(N); % interaction matrix
Ki=1*ones(N,1); % death rate
t0=0; T=300; % t0 is start time, and T is the end time
b=[4, .95, 1.05]; % growth rate

alpha=[.6,.6,1]; % order of derivatives
% alpha=[1,1,1]; % order of derivatives

h=0.01; % step size for computing
for i=1:2*NXB0
    if i>NXB0
    x0=[XB0(i-NXB0),.1,1]; % initial conditions
    else
    x0=[XB0(i),.3,.1]; % initial conditions
    end
[t, x] = FDE_PI12_PC(alpha,@fun,t0,T,x0',h);
RelX=x./(ones(N,1)*sum(x)); % Relative abundances

XB(i,:,:)=RelX;
end
%% panel b

f=figure;
f.Position = [0 0 1800 500];
f.Renderer='painters';
% defining blindfriendly colors (red and green)
PcR= [0.92,0.27,0.18];
PcG= [0.18,0.40,0.14];
   
%        ------------------------------------------------------------------    
    % Plot the  abundance of species
    hold on
    Xb(:,:)=XB(:,1,1:20:end);
    Xr(:,:)=XB(:,2,1:20:end);
    Xg(:,:)=XB(:,3,1:20:end);
    tt=t(:,1:20:end);
subplot(1,3,1)    
for j=1:i
    if j==2 || j==8
%             p=semilogy(tt,Xb(j,:),'--b');
                       
    else
%     p=semilogy(tt,Xb(j,:),'b');p.Color(4)=.3;
    p=plot(tt,Xb(j,:),'b');p.Color(4)=.3;
    end
    hold on
    set(p,'LineWidth',2)
end
 plot(tt,Xb(2,:),'--b','LineWidth',2);
 p=plot(tt,Xb(8,:),':b','LineWidth',3);
    hold on
%     p3=semilogy(tt(1:100:end),Xb(7,1:100:end),'ob','MarkerSize', 10);
%     p4=semilogy(tt(1:500:end),Xb(6,1:500:end),'sqb','MarkerSize', 10);
    ylabel('Abundance', 'FontSize',28)
%     
%     yticks([0.01 0.1 1])
% yticklabels({'-2','-1','0'})
% Settings for the general plot

xlabel('Time')
axis([0,T,.001,1])
set(gca,'FontSize',20)
% set(gca, 'YScale', 'log')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,2)

hold on
    for j=1:i
        if j==2 || j==8
%             p1=semilogy(tt,Xr(j,:),'--','color',PcR);
           
        else
                p1=plot(tt,Xr(j,:),'color',PcR);p1.Color(4)=.3;
            %     p1=semilogy(tt,Xr(j,:),'color',PcR);p1.Color(4)=.3;
        end
    set(p1,'LineWidth',2)
    end
     plot(tt,Xr(2,:),'--','color',PcR,'LineWidth',2);
          p1=plot(tt,Xr(8,:),':','color',PcR,'LineWidth',3);
%     semilogy(tt(1:100:end),Xr(7,1:100:end),'o','color',PcR,'MarkerSize', 10);
%     semilogy(tt(1:100:end),Xr(6,1:100:end),'sq','MarkerSize', 10,'LineWidth',2);
    
%     
%     yticks([0.01 0.1 1])
% yticklabels({'-2','-1','0'})
xlabel('Time')
axis([0,T,.001,1])
set(gca,'FontSize',20)
box on
% set(gca, 'YScale', 'log')
%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,3)

hold on
    for j=1:i
        if j==2 || j==8

%             p2=semilogy(tt,Xg(j,:),'--','color',PcG);
        else
%     p2=semilogy(tt,Xg(j,:),'color',PcG);p2.Color(4)=.3;
p2=plot(tt,Xg(j,:),'color',PcG);p2.Color(4)=.3;
        end
    set(p2,'LineWidth',2)
    end
        plot(tt,Xg(2,:),'--','color',PcG,'LineWidth',2);
                p2=plot(tt,Xg(8,:),':','color',PcG,'LineWidth',3);
%     semilogy(tt(1:100:end),Xg(7,1:100:end),'om','MarkerSize', 10);
%     semilogy(tt(1:100:end),Xg(6,1:100:end),'sqk','MarkerSize', 10,'LineWidth',2);
    
    
%     yticks([0.01 0.1 1])
% yticklabels({'-2','-1','0'})
xlabel('Time')
axis([0,T,.001,1])
set(gca,'FontSize',20)
box on
% set(gca, 'YScale', 'log')
%%
alpha=[1,1,1]; % order of derivatives
T=300;
h=0.01; % step size for computing
for i=1:2*NXB0
    if i>NXB0
    x0=[XB0(i-NXB0),.1,1]; % initial conditions
    else
    x0=[XB0(i),.3,.1]; % initial conditions
    end
[t, x] = FDE_PI12_PC(alpha,@fun,t0,T,x0',h);
RelX1=x./(ones(N,1)*sum(x)); % Relative abundances

XB1(i,:,:)=RelX1;
end
%% Panel a 

f=figure;
f.Position = [0 0 1800 500];
f.Renderer='painters';
% defining blindfriendly colors (red and green)
PcR= [0.92,0.27,0.18];
PcG= [0.18,0.40,0.14];
   
%        ------------------------------------------------------------------    
    % Plot the  abundance of species
    hold on
    Xb1(:,:)=XB1(:,1,1:20:end);
    Xr1(:,:)=XB1(:,2,1:20:end);
    Xg1(:,:)=XB1(:,3,1:20:end);
    tt=t(:,1:20:end);
subplot(1,3,1)    
for j=1:i
    if j==2 || j==8
%             p=semilogy(tt,Xb(j,:),'--b');
                       
    else
%     p=semilogy(tt,Xb(j,:),'b');p.Color(4)=.3;
    p=plot(tt,Xb1(j,:),'b');p.Color(4)=.3;
    end
    hold on
    set(p,'LineWidth',2)
end
 plot(tt,Xb1(2,:),'--b','LineWidth',2);
 p=plot(tt,Xb1(8,:),':b','LineWidth',3);
    hold on
%     p3=semilogy(tt(1:100:end),Xb(7,1:100:end),'ob','MarkerSize', 10);
%     p4=semilogy(tt(1:500:end),Xb(6,1:500:end),'sqb','MarkerSize', 10);
    ylabel('Abundance', 'FontSize',28)
%     
%     yticks([0.01 0.1 1])
% yticklabels({'-2','-1','0'})
% Settings for the general plot

xlabel('Time')
axis([0,T,.001,1])
set(gca,'FontSize',20)
% set(gca, 'YScale', 'log')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,2)

hold on
    for j=1:i
        if j==2 || j==8
%             p1=semilogy(tt,Xr(j,:),'--','color',PcR);
           
        else
                p1=plot(tt,Xr1(j,:),'color',PcR);p1.Color(4)=.3;
            %     p1=semilogy(tt,Xr(j,:),'color',PcR);p1.Color(4)=.3;
        end
    set(p1,'LineWidth',2)
    end
     plot(tt,Xr1(2,:),'--','color',PcR,'LineWidth',2);
          p1=plot(tt,Xr1(8,:),':','color',PcR,'LineWidth',3);
%     semilogy(tt(1:100:end),Xr(7,1:100:end),'o','color',PcR,'MarkerSize', 10);
%     semilogy(tt(1:100:end),Xr(6,1:100:end),'sq','MarkerSize', 10,'LineWidth',2);
    
%     
%     yticks([0.01 0.1 1])
% yticklabels({'-2','-1','0'})
xlabel('Time')
axis([0,T,.001,1])
set(gca,'FontSize',20)
box on
% set(gca, 'YScale', 'log')
%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,3)

hold on
    for j=1:i
        if j==2 || j==8

%             p2=semilogy(tt,Xg(j,:),'--','color',PcG);
        else
%     p2=semilogy(tt,Xg(j,:),'color',PcG);p2.Color(4)=.3;
p2=plot(tt,Xg1(j,:),'color',PcG);p2.Color(4)=.3;
        end
    set(p2,'LineWidth',2)
    end
        plot(tt,Xg1(2,:),'--','color',PcG,'LineWidth',2);
                p2=plot(tt,Xg1(8,:),':','color',PcG,'LineWidth',3);
%     semilogy(tt(1:100:end),Xg(7,1:100:end),'om','MarkerSize', 10);
%     semilogy(tt(1:100:end),Xg(6,1:100:end),'sqk','MarkerSize', 10,'LineWidth',2);
    
    
%     yticks([0.01 0.1 1])
% yticklabels({'-2','-1','0'})
xlabel('Time')
axis([0,T,.001,1])
set(gca,'FontSize',20)
box on
% set(gca, 'YScale', 'log')
