%% Main_ThreeSpecies
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
global n N Ki b Kij T mu x0
%% Coefficients and Conditions

mu=[1,1,1]; % [mu_B,mu_R,mu_G] Order of derivatives,  0<mu(i)=<1

n=2; % Hill coefficient

N=3; % number of species

Kij=0.1*ones(N); % interaction matrix

Ki=1*ones(N,1); % death rate

T=150; %  final time

Perturbation='OUP'; % Possible usages 'False', 'Pulse1', 'Pulse2', 'Pulse3', 'Pulse4', 'Periodic', 'OUP', and 'OUP_new'

b=[1, .95, 1.05]; % growth rates for cases: False, Pulse, and Periodic

x0=[.3,.1,.2]'; % initial conditions

[t,x,B]=method1(Perturbation);


%% plotting

RelX=x./(ones(N,1)*sum(x)); % Relative abundances

f=figure;
f.Renderer='painters';

% defining blindfriendly colors (red and green)
PcR= [0.92,0.27,0.18];
PcG= [0.18,0.40,0.14];

switch Perturbation    
    case 'False'        

%        ------------------------------------------------------------------    
    % Plot the abundance of species
    p=plot(t,x(1,:),'b',t,x(2,:),'r',t,x(3,:),'g');
    ylabel('Species abundance')

%        ------------------------------------------------------------------    

    case 'Pulse1'
        
    %%Plot of growth rate with pulse        
        
    tt=0:.1:T;Nt=length(tt); % time simulating with 0.1 step size
        
    % Simulate growth rates with pulse
    BB=zeros(1,Nt);RR=zeros(1,Nt);GG=zeros(1,Nt);
    for i=1:Nt
        if tt(i)>20 && tt(i)<60    
        BB(i)=.5; GG(i)=2;    
        else
        BB(i)=1; GG(i)=1.05;
        end
        RR(i)=b(2);
    end
    
    % Highlight background as perturbation
    vb = [20 .5; 60 .5; 60 2.2; 20 2.2];
    f = [1 2 3 4];
    h11=patch('Faces',f,'Vertices',vb,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    set(gca,'YScale','log')
    hold on
    
    % plot (log) of growth rates    
    pb=semilogy(tt,BB,'b',tt,RR,'r',tt,GG,'g');
    
    % Settings for plot
    set(pb,'LineWidth',4)
    ylabel('Growth Rate (log)')
    xlabel('Time')
    legend('Pulse: b_{B}=0.5, b_{G}=2','b_B\approx1','b_R=0.95', 'b_G\approx1.05')
    set(gca,'FontSize',23, 'FontWeight', 'bold')
    
    %%plotting relative abundance of species
    figure
    v = [20 0; 60 0; 60 1; 20 1];
    h1=patch('Faces',f,'Vertices',v,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
     
    hold on
    
    % Plot the relative abundances
    p=plot(t,RelX(1,:),'b',t,RelX(2,:),'r',t,RelX(3,:),'g');
    
    % Remove highlighted bars from the legend
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
       ylabel('Relative abundance')
       
%        ------------------------------------------------------------------    

    case 'Pulse2'
        
    %%Plot of growth rate with pulse        
        
    tt=0:.1:T;Nt=length(tt); % time simulating with 0.1 step size
        
    % Simulate growth rates with pulse
    BB=zeros(1,Nt);RR=zeros(1,Nt);GG=zeros(1,Nt);
    for i=1:Nt
        if tt(i)>20 && tt(i)<60    
        BB(i)=.5; GG(i)=2.2;    
        else
        BB(i)=1; GG(i)=1.05;
        end
        RR(i)=b(2);
    end
    
    % Highlight background as perturbation
    vb = [20 .5; 60 .5; 60 2.2; 20 2.2];
    f = [1 2 3 4];
    h11=patch('Faces',f,'Vertices',vb,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    set(gca,'YScale','log')
    hold on
    
    % plot (log) of growth rates    
    pb=semilogy(tt,BB,'b',tt,RR,'r',tt,GG,'g');
    
    % Settings for plot
    set(pb,'LineWidth',4)
    ylabel('Growth Rate (log)')
    xlabel('Time')
    legend('Pulse: b_{B}=0.5, b_{G}=2','b_B\approx1','b_R=0.95', 'b_G\approx1.05')
    set(gca,'FontSize',23, 'FontWeight', 'bold')
    
    %%plotting relative abundance of species
    figure
    v = [20 0; 60 0; 60 1; 20 1];
    h1=patch('Faces',f,'Vertices',v,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
     
    hold on
    
    % Plot the relative abundances
    p=plot(t,RelX(1,:),'b',t,RelX(2,:),'r',t,RelX(3,:),'g');
    
    % Remove highlighted bars from the legend
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
       ylabel('Relative abundance')
       
%        ------------------------------------------------------------------    

    case 'Pulse3'
        
    %%Plot of growth rate with pulse        
        
    tt=0:.1:T;Nt=length(tt); % time simulating with 0.1 step size
        
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
    h11=patch('Faces',f,'Vertices',vb,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
    v1b = [200 0.2; 330 0.2; 330 4.5; 200 4.5];
    h22=patch('Faces',f,'Vertices',v1b,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.36);
    set(gca,'YScale','log')
    
    % plot (log) of growth rates    
    pb=semilogy(tt,BB,'b',tt,RR,'r',tt,GG,'g');
    
    % Settings for plot
    set(pb,'LineWidth',4)
    ylabel('Growth Rate (log)')
    xlabel('Time')
    legend('Perturbation1: b_{B}=0.2','Perturbation2: b_{B}=4.5','b_B\approx1','b_R=0.95', 'b_G=1.05')
    set(gca,'FontSize',23, 'FontWeight', 'bold')
    
    %%plotting relative abundance of species
    figure
    v = [60 0; 100 0; 100 1; 60 1];
    v1 = [200 0; 330 0; 330 1; 200 1];
    h1=patch('Faces',f,'Vertices',v,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
    h2=patch('Faces',f,'Vertices',v1,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.36);
    
    % Plot the relative abundances
    p=plot(t,RelX(1,:),'b',t,RelX(2,:),'r',t,RelX(3,:),'g');
    
    % Remove highlighted bars from the legend
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
       ylabel('Relative abundance')
%        ------------------------------------------------------------------    

    case 'Pulse4'
        
    %%Plot of growth rate with pulse        
        
    tt=0:.1:T;Nt=length(tt); % time simulating with 0.1 step size
        
    % Simulate growth rates with pulse
    BB=zeros(1,Nt);RR=zeros(1,Nt);GG=zeros(1,Nt);
    for i=1:Nt
        if tt(i)>60 && tt(i)<100    
        BB(i)=.2;    
        elseif tt(i)>400 && tt(i)<530    
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
    h11=patch('Faces',f,'Vertices',vb,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
    v1b = [400 0.2; 530 0.2; 530 4.5; 400 4.5];
    h22=patch('Faces',f,'Vertices',v1b,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.36);
    set(gca,'YScale','log')
    
    % plot (log) of growth rates    
    pb=semilogy(tt,BB,'b',tt,RR,'r',tt,GG,'g');
    
    % Settings for plot
    set(pb,'LineWidth',4)
    ylabel('Growth Rate (log)')
    xlabel('Time')
    legend('Perturbation1: b_{B}=0.2','Perturbation2: b_{B}=4.5','b_B\approx1','b_R=0.95', 'b_G=1.05')
    set(gca,'FontSize',23, 'FontWeight', 'bold')
    
    %%plotting relative abundance of species
    figure
    v = [60 0; 100 0; 100 1; 60 1];
    v1 = [200 0; 330 0; 330 1; 200 1];
    h1=patch('Faces',f,'Vertices',v,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
    h2=patch('Faces',f,'Vertices',v1,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.36);
    
    % Plot the relative abundances
    p=plot(t,RelX(1,:),'b',t,RelX(2,:),'r',t,RelX(3,:),'g');
    
    % Remove highlighted bars from the legend
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
       ylabel('Relative abundance')

%        ------------------------------------------------------------------
       
    case 'Periodic'
 
    %%plot of growth rate with periodic perturbation
    
    tt=0:.1:T;Nt=length(tt); % time simulating with 0.1 step size 
    
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
    f=1:4*m;f=reshape(f,4,m);f=f';
    h11=patch('Faces',f,'Vertices',vb,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
    h22=patch('Faces',f,'Vertices',v1b,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.36);   
    set(gca,'YScale','log')
    
    % plot (log) of growth rates
    pb=semilogy(tt,BB,'b',tt,RR,'r',tt,GG,'g');
    
    % Settings for plot
    set(pb,'LineWidth',4)
    ylabel('Growth Rate (log)')
    xlabel('Time')
    legend('Perturbation1: b_{B}=0.2','Perturbation2: b_{B}=4.5','b_B\approx1','b_R=0.95', 'b_G=1.05')
    set(gca,'FontSize',23, 'FontWeight', 'bold')
    
    %%plotting relative abundance of species
    figure        
    
    % Highlight background as perturbation
    m=ceil(mod(T/(mm*4),mm));
    v=zeros(4*m,2);
    v1=zeros(4*m,2);
    for i=1:m    
        v(4*i-3:i*4,:)=[mm*(4*i-3) 0; mm*(4*i-2) 0; mm*(4*i-2) 1; mm*(4*i-3) 1];
        v1(4*i-3:i*4,:)=[mm*(4*i-1) 0; mm*(4*i) 0; mm*(4*i) 1; mm*(4*i-1) 1];
    end
    h1=patch('Faces',f,'Vertices',v,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
    h2=patch('Faces',f,'Vertices',v1,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.36);
    
    % Plot the relative abundances
    p=plot(t,RelX(1,:),'b',t,RelX(2,:),'r',t,RelX(3,:),'g');
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
       ylabel('Relative abundance')

%        ------------------------------------------------------------------


    otherwise
        
    %%plot of growth rate with stochastic perturbation
    
    tt=1:T; % time simulating with 1 step size  
    
    % Plot the simulated growth rates
    pb=plot(tt,B(1:T,1),'b',tt,B(1:T,2),'r',tt,B(1:T,3),'g');
    
    % Settings for plot
    set(pb,'LineWidth',4)
    ylabel('Growth Rates: OUP')
    xlabel('Time')
    legend('b_B','b_R', 'b_G')
    set(gca,'FontSize',23, 'FontWeight', 'bold')
    
    %%plotting relative abundance of species
       figure
            p=plot(t,RelX(1,:),'b',t,RelX(2,:),'r',t,RelX(3,:),'g');
       ylabel('Relative abundance')
end

% Settings for the general plot
set(p,'LineWidth',5)
set(gca,'FontSize',23, 'FontWeight', 'bold')
xlabel('Time')

% % Generate legend for showing memory of each species
legend(['X_{B}/X_{total}, Memory_B=',num2str(1-mu(1))],...
    ['X_{R}/X_{total}, Memory_R=',num2str(1-mu(2))],...
    ['X_{G}/X_{total}, Memory_G=',num2str(1-mu(3))],'Location', 'Best');



