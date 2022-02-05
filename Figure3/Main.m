%% Figure 3
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
%        b - Growth rates (including pulse and Periodic perturbations), e.g. b=[1, .95, 1.05];
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
global n N Ki b Kij x0
%% Inputs
% Coefficients and Conditions

% [mu_B,mu_R,mu_G] Order of derivatives,  0<mu(i)=<1
mu1=[1,1,1]; % for panel b and d, row 1 (no memoey) 
mu2=0.96*[1,1,1]; % for panel b and d, row 2 (with memoey)
mu3=0.9*[1,1,1]; % for panel b and d, row 3 (with memoey)

n=2; % Hill coefficient

N=3; % number of species

Kij=0.1*ones(N); % interaction matrix

Ki=1*ones(N,1); % death rate

T1=450; %  final time
T2=720; %  final time

b=[1, .95, 1.05]; % growth rates for cases: False, Pulse, and Periodic

x0=[.99;.01;.01]; % initial conditions

t0=0; % initial time
h=0.01; % step size for computing

Fun1=@fun13; % ODE model including a pulse for panel b 
Fun2=@fun2;  % ODE model including a pulse for panel d

% solver for fractional differential equation
[t1, x1] = FDE_PI12_PC(mu1,Fun1,t0,T1,x0,h); % results for panel b row 1
[t2, x2] = FDE_PI12_PC(mu1,Fun2,t0,T2,x0,h); % results for panel d row 1
[~, x12] = FDE_PI12_PC(mu2,Fun1,t0,T1,x0,h); % results for panel b row 2
[~, x22] = FDE_PI12_PC(mu2,Fun2,t0,T2,x0,h); % results for panel d row 2
[~, x13] = FDE_PI12_PC(mu3,Fun1,t0,T1,x0,h); % results for panel b row 3
[~, x23] = FDE_PI12_PC(mu3,Fun2,t0,T2,x0,h); % results for panel d row 3
%% plotting

RelX1=x1./(ones(N,1)*sum(x1)); % Relative abundances
RelX2=x2./(ones(N,1)*sum(x2));
RelX12=x12./(ones(N,1)*sum(x12));
RelX22=x22./(ones(N,1)*sum(x22));
RelX13=x13./(ones(N,1)*sum(x13));
RelX23=x23./(ones(N,1)*sum(x23));

% defining blindfriendly colors (red and green)
PcR= [0.92,0.27,0.18];
PcG= [0.18,0.40,0.14];
            
%% Plot of growth rate with pulse (panel a)
  f1=figure;
f1.Renderer='painters';

    tt=0:.1:T1;Nt=length(tt); % time simulating with 0.1 step size
        
    % Simulate growth rates with pulse
    BB=zeros(1,Nt);RR=zeros(1,Nt);GG=zeros(1,Nt);
    for i=1:Nt
        if tt(i)>60 && tt(i)<100    
        BB(i)=.2;    
        elseif tt(i)>200 && tt(i)<330    
        BB(i)=4.5;
        else
        BB(i)=1;
        end
        RR(i)=b(2);
        GG(i)=b(3);
    end
    
    % Highlight background as perturbation
    vb = [60 .2; 100 .2; 100 4.5; 60 4.5];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices',vb,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
    v1b = [200 0.2; 330 0.2; 330 4.5; 200 4.5];
    patch('Faces',f,'Vertices',v1b,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.36);
    set(gca,'YScale','log')
    
    semilogy(tt,BB,'Color',[0,0,1],'LineWidth',4);
    hold on
    semilogy(tt,RR,'Color',PcR,'LineWidth',4);
    semilogy(tt,GG,'Color',PcG,'LineWidth',4);
%     Settings for plot
axis tight
    
    % Settings for plot
    ylabel('Log growth rate')
    xlabel('Time')
    legend('Perturbation1: b_{B}=0.2','Perturbation2: b_{B}=4.5','b_B\approx1','b_R=0.95', 'b_G=1.05')
    set(gca,'FontSize',23, 'FontWeight', 'bold')

%%  Plot of growth rate with 'Periodic' perturbations (panel c)
 f12=figure;
 f12.Renderer='painters';
    %%plot of growth rate with periodic perturbation
    
    tt=0:.1:T2;Nt=length(tt); % time simulating with 0.1 step size 
    
    % Simulate growth rates with periodic perturbation
    mm=20;
    BB=zeros(1,Nt);RR=zeros(1,Nt);GG=zeros(1,Nt);
    for i=1:Nt
        m=ceil(mod(tt(i)/(mm*4),mm));
        if tt(i)>mm*(4*m-3) && tt(i)<mm*(4*m-2)    
            BB(i)=.2;    
        elseif tt(i)>mm*(4*m-1) && tt(i)<mm*(4*m)
            BB(i)=4.5;
        else
            BB(i)=1;
        end
        RR(i)=b(2);
        GG(i)=b(3);
    end
    
    % Highlight background as perturbation
    vb=zeros(4*m,2);
    v1b=zeros(4*m,2);
    for i=1:m    
        vb(4*i-3:i*4,:)=[mm*(4*i-3) .2; mm*(4*i-2) .2; mm*(4*i-2) 4.5; mm*(4*i-3) 4.5];
        v1b(4*i-3:i*4,:)=[mm*(4*i-1) .2; mm*(4*i) .2; mm*(4*i) 4.5; mm*(4*i-1) 4.5];
    end
    ff=1:4*m;ff=reshape(ff,4,m);ff=ff';
    h11=patch('Faces',ff,'Vertices',vb,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
    h22=patch('Faces',ff,'Vertices',v1b,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.36);   
    set(gca,'YScale','log')
    
    % plot (log) of growth rates
     
    semilogy(tt,BB,'Color',[0,0,1],'LineWidth',4);
    hold on
    semilogy(tt,RR,'Color',PcR,'LineWidth',4);
    semilogy(tt,GG,'Color',PcG,'LineWidth',4);
    axis tight

    % Settings for plot
    ylabel('Log growth rates')
    xlabel('Time')
    legend('Perturbation1: b_{B}=0.2','Perturbation2: b_{B}=4.5','b_B\approx1','b_R=0.95', 'b_G=1.05')
    set(gca,'FontSize',23, 'FontWeight', 'bold')

%% plotting relative abundance of species for panel b (pulse) row 1
 
  f13=figure;
  f13.Renderer='painters';

    v = [60 0; 100 0; 100 1; 60 1];
    v1 = [200 0; 330 0; 330 1; 200 1];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices',v,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
    patch('Faces',f,'Vertices',v1,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.36);
    
    % Plot the relative abundances
    semilogy(t1,RelX1(1,:),'Color',[0,0,1],'LineWidth',4);
    hold on
    semilogy(t1,RelX1(2,:),'Color',PcR,'LineWidth',4);
    semilogy(t1,RelX1(3,:),'Color',PcG,'LineWidth',4);
    
    % Remove highlighted bars from the legend
%     set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
       ylabel('Relative abundance')    
     xlabel('Time')
       set(gca,'FontSize',17, 'FontWeight', 'bold')
%% plotting relative abundance of species for panel b (pulse) row 2
 
  f132=figure;
  f132.Renderer='painters';

    v = [60 0; 100 0; 100 1; 60 1];
    v1 = [200 0; 330 0; 330 1; 200 1];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices',v,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
    patch('Faces',f,'Vertices',v1,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.36);
    
    % Plot the relative abundances
    semilogy(t1,RelX12(1,:),'Color',[0,0,1],'LineWidth',4);
    hold on
    semilogy(t1,RelX12(2,:),'Color',PcR,'LineWidth',4);
    semilogy(t1,RelX12(3,:),'Color',PcG,'LineWidth',4);
    
    % Remove highlighted bars from the legend
%     set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
       ylabel('Relative abundance')    
     xlabel('Time')
       set(gca,'FontSize',17, 'FontWeight', 'bold')
%% plotting relative abundance of species for panel b (pulse) row 3
 
  f133=figure;
  f133.Renderer='painters';

    v = [60 0; 100 0; 100 1; 60 1];
    v1 = [200 0; 330 0; 330 1; 200 1];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices',v,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
    patch('Faces',f,'Vertices',v1,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.36);
    
    % Plot the relative abundances
    semilogy(t1,RelX13(1,:),'Color',[0,0,1],'LineWidth',4);
    hold on
    semilogy(t1,RelX13(2,:),'Color',PcR,'LineWidth',4);
    semilogy(t1,RelX13(3,:),'Color',PcG,'LineWidth',4);
    
    % Remove highlighted bars from the legend
%     set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
       ylabel('Relative abundance')    
     xlabel('Time')
       set(gca,'FontSize',17, 'FontWeight', 'bold')
    
%% plotting relative abundance of species for panel d (periodic) row 1
    f10=figure;
f10.Renderer='painters';        
    
    % Highlight background as perturbation
    m=ceil(mod(T2/(mm*4),mm));
    v=zeros(4*m,2);
    v1=zeros(4*m,2);
    for i=1:m    
        v(4*i-3:i*4,:)=[mm*(4*i-3) 0; mm*(4*i-2) 0; mm*(4*i-2) 1; mm*(4*i-3) 1];
        v1(4*i-3:i*4,:)=[mm*(4*i-1) 0; mm*(4*i) 0; mm*(4*i) 1; mm*(4*i-1) 1];
    end
    h1=patch('Faces',ff,'Vertices',v,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
    h2=patch('Faces',ff,'Vertices',v1,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.36);
    
    % Plot the relative abundances
    semilogy(t2,RelX2(1,:),'Color',[0,0,1],'LineWidth',4);
    hold on
    semilogy(t2,RelX2(2,:),'Color',PcR,'LineWidth',4);
    semilogy(t2,RelX2(3,:),'Color',PcG,'LineWidth',4);
%     
%     
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
       ylabel('Relative abundance')

% Settings for the general plot
set(gca,'FontSize',23, 'FontWeight', 'bold')
xlabel('Time')

axis tight 

% % Generate legend for showing memory of each species
legend(['X_{B}/X_{total}, Memory_B=',num2str(1-mu1(1))],...
    ['X_{R}/X_{total}, Memory_R=',num2str(1-mu1(2))],...
    ['X_{G}/X_{total}, Memory_G=',num2str(1-mu1(3))],'Location', 'Best');

%% plotting relative abundance of species for panel d (periodic) row 2
    f14=figure;
f14.Renderer='painters';        
    
    % Highlight background as perturbation
    m=ceil(mod(T2/(mm*4),mm));
    v=zeros(4*m,2);
    v1=zeros(4*m,2);
    for i=1:m    
        v(4*i-3:i*4,:)=[mm*(4*i-3) 0; mm*(4*i-2) 0; mm*(4*i-2) 1; mm*(4*i-3) 1];
        v1(4*i-3:i*4,:)=[mm*(4*i-1) 0; mm*(4*i) 0; mm*(4*i) 1; mm*(4*i-1) 1];
    end
    h1=patch('Faces',ff,'Vertices',v,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
    h2=patch('Faces',ff,'Vertices',v1,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.36);
    
    % Plot the relative abundances
    semilogy(t2,RelX22(1,:),'Color',[0,0,1],'LineWidth',4);
    hold on
    semilogy(t2,RelX22(2,:),'Color',PcR,'LineWidth',4);
    semilogy(t2,RelX22(3,:),'Color',PcG,'LineWidth',4);
%     
%     
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
       ylabel('Relative abundance')

% Settings for the general plot
set(gca,'FontSize',23, 'FontWeight', 'bold')
xlabel('Time')

axis tight 

% % Generate legend for showing memory of each species
legend(['X_{B}/X_{total}, Memory_B=',num2str(1-mu2(1))],...
    ['X_{R}/X_{total}, Memory_R=',num2str(1-mu2(2))],...
    ['X_{G}/X_{total}, Memory_G=',num2str(1-mu2(3))],'Location', 'Best');
%% plotting relative abundance of species for panel d (periodic) row 3
    f15=figure;
f15.Renderer='painters';        
    
    % Highlight background as perturbation
    m=ceil(mod(T2/(mm*4),mm));
    v=zeros(4*m,2);
    v1=zeros(4*m,2);
    for i=1:m    
        v(4*i-3:i*4,:)=[mm*(4*i-3) 0; mm*(4*i-2) 0; mm*(4*i-2) 1; mm*(4*i-3) 1];
        v1(4*i-3:i*4,:)=[mm*(4*i-1) 0; mm*(4*i) 0; mm*(4*i) 1; mm*(4*i-1) 1];
    end
    h1=patch('Faces',ff,'Vertices',v,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
    h2=patch('Faces',ff,'Vertices',v1,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.36);
    
    % Plot the relative abundances
    semilogy(t2,RelX23(1,:),'Color',[0,0,1],'LineWidth',4);
    hold on
    semilogy(t2,RelX23(2,:),'Color',PcR,'LineWidth',4);
    semilogy(t2,RelX23(3,:),'Color',PcG,'LineWidth',4);
%     
%     
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
       ylabel('Relative abundance')

% Settings for the general plot
set(gca,'FontSize',23, 'FontWeight', 'bold')
xlabel('Time')

axis tight 

% % Generate legend for showing memory of each species
legend(['X_{B}/X_{total}, Memory_B=',num2str(1-mu3(1))],...
    ['X_{R}/X_{total}, Memory_R=',num2str(1-mu3(2))],...
    ['X_{G}/X_{total}, Memory_G=',num2str(1-mu3(3))],'Location', 'Best');
