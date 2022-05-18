%% Figure 7 panel C
clear 
clc

global A mu

order1=[1,1];

mu=[0.599 0.626];

t0=0;
T=300;
h=.1;
F=@fun;
JF=@Jfun;

A=[-0.9059 -0.9377;-0.972 -0.9597];

%% fix points
xx1=[(A(1,2)*mu(2)-A(2,2)*mu(1))/(A(1,1)*A(2,2)-A(1,2)*A(2,1));...
    (A(2,1)*mu(1)-A(1,1)*mu(2))/(A(1,1)*A(2,2)-A(1,2)*A(2,1))];
xx3=[1e-3;-mu(2)/A(2,2)];
X0=xx3; % initial conditions

p=.2; %perturb
[t,X]=FDE_PI2_IM(order1,F,JF,t0,T,X0,h,p);

%%

        RelX=X./(ones(2,1)*sum(X));
    %%Plot of growth rate with pulse        
        
    tt=0:.1:T;Nt=length(tt); % time simulating with 0.1 step size
        
    % Simulate growth rates with pulse
    BB=zeros(1,Nt);RR=zeros(1,Nt);
    for i=1:Nt
        if tt(i)>50 && tt(i)<80  
        BB(i)=mu(1)+p;     
        else
        BB(i)=mu(1); 
        end
        RR(i)=mu(2);
    end
    
    % Highlight background as perturbation
    vb = [50 0.5; 80 0.5; 80 .9; 50 .9];
    f = [1 2 3 4];
    h11=patch('Faces',f,'Vertices',vb,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
%     set(gca,'YScale','log')
    hold on
    
    % plot (log) of growth rates    
    pb=plot(tt,BB,'b',tt,RR,'r');
    
    % Settings for plot
    set(pb,'LineWidth',4)
    ylabel('Growth Rate')
    xlabel('Time')   
annotation('textarrow',[0.48 0.38],[0.67 0.72],'String','perturbation + b_{BT}', 'FontSize',19)
    set(gca,'FontSize',23)
    
    axis tight
    
    %%plotting relative abundance of species
    figure
    v = [50 0; 80 0; 80 1; 50 1];
    h1=patch('Faces',f,'Vertices',v,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
     
    hold on
    
    % Plot the relative abundances
    p=plot(t,RelX(1,:),'b',t,RelX(2,:),'r');
    set(p,'LineWidth',4)
    % Remove highlighted bars from the legend
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
       ylabel('Abundance')
       
        set(gca,'FontSize',23)
    xlabel('Time') 
       axis tight
%        set(gcf,'renderer','Painters')
