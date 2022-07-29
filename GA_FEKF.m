%Guel-Cortez et al 2022. 
%Parameter Estimation of Fractional-Order Systems via Evolutionary Algorithms and the Extended Fractional Kalman Filter
%General script for fractional EKF with GA in a mass spring damper system

close all;
clearvars;
clc

timerVal=tic;
rng('default');
rng(1); %1 error 1.0893

L = 100;
Tspan=25;
T=.05;
t=0:T:Tspan;
ks=1.5;
b=0.8;
N=[.8;.8];
R=1e-4;
Q=1e-6*eye(2,2);
u=5*sin(t);
[t,y]=FOsystem(ks,b,N,Q,R,T,t,L,u);
Q=1e-6*eye(5,5); Q(3,3)=0; Q(4,4)=0;Q(5,5)=0;
realp=[ks,b,N(1)];

options = optimoptions('ga','PlotFcn',"gaplotbestf",'UseParallel', true,'MaxStallGenerations',500,'MaxGenerations',1000);
[x,fval] = ga(@(gains) FOKFilter(t,u,y,gains,Q,R,T,L),3,[],[],[],[],[-1;-1;0],[2;2;1],[],options);

[RMSE,te,ye,x_e]=FOKFilter(t,u,y,x,Q,R,T,L);
figure 
set(gcf,'color','w');
plot(t,y,te,ye)
xlabel('t (s)','Interpreter','Latex','FontSize', 12)
ylabel('$y(t)/x_e(t)$','Interpreter','Latex','FontSize', 12)
legend('Measurement','Response with estimated parameters')
% 
 figure
   set(gcf,'color','w');
subplot(3,2,[1 2])
plot(t,y(1,:),t,x_e(1,:))
xlabel('t (s)','Interpreter','Latex','FontSize', 12)
ylabel('$y(t)/x_e(t)$','Interpreter','Latex','FontSize', 12)
legend('Measurement','Estimated')
subplot(3,2,3)
plot(t,x_e(2,:))
xlabel('t (s)','Interpreter','Latex','FontSize', 12)
ylabel('$\dot{x}_e(t)$','Interpreter','Latex','FontSize', 12)
subplot(3,2,4)
plot(t,x_e(3,:))
xlabel('t (s)','Interpreter','Latex','FontSize', 12)
ylabel('$k_e(t)$','Interpreter','Latex','FontSize', 12)
subplot(3,2,5)
plot(t,x_e(4,:))
xlabel('t (s)','Interpreter','Latex','FontSize', 12)
ylabel('$b_e(t)$','Interpreter','Latex','FontSize', 12)
subplot(3,2,6)
plot(t,x_e(5,:))
%ylim([N(1)-0.1,N(1)+0.1])
xlabel('t (s)','Interpreter','Latex','FontSize', 12)
ylabel('$\alpha_e(t)$','Interpreter','Latex','FontSize', 12)

estimatedp=[x_e(3:5,end)]
norm(realp-estimatedp,2)
Elapsetime=toc(timerVal);
