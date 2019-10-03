clear;clc
addpath Tool/
syms A(t) Ms(t) Me(t) Mi(t) Hs(t) He(t) Hi(t) Hr(t) Hit(t)
syms delta C gamma_m mu_a f beta_m mu_m theta_m mu_h beta_h theta_h gamma_h M H
M=Ms+Me+Mi;
H=Hs+He+Hi+Hr;
ode1 = diff(A) == delta*(1-A/C)*M - (gamma_m+mu_a)*A;
ode2 = diff(Ms) == f*gamma_m*A - beta_m*Hi/H*Ms - mu_m*Ms;
ode3 = diff(Me) ==  beta_m*Hi/H*Ms - (theta_m+mu_m)*Me;
ode4 = diff(Mi) == theta_m*Me - mu_m*Mi;
ode5 = diff(Hs) == mu_h*H - beta_h*Mi/M*Hs - mu_h*Hs;
ode6 = diff(He) == beta_h*Mi/M*Hs - (mu_h+theta_h)*He;
ode7 = diff(Hi) == theta_h*He - (gamma_h+mu_h)*Hi;
ode8 = diff(Hr) == gamma_h*Hi-mu_h*Hr;
ode9 = diff(Hit) == theta_h*He;
odes=[ode1; ode2; ode3; ode4; ode5; ode6; ode7; ode8; ode9];
vars=[Hit Hi A Ms Me Mi Hs He Hr];

% init=[9000 1199950 40 10 321710 18 6 81501];
% pars=[10000 1.3 0.15 0.6 65 0.5 1 0.9 0.13 0.0004 0.12 0.7 0.6];
% range=[5755, 17265; 0 8e6; 0 100; 0 100; 244402 321734; 18 72; ...
%     6 24; 81405 158809; 6400 95000; 1 1.5; 0 4; 0 4; 55 165; 0.42 0.55; 0.5 1.75; 0.875 1.4; 0.007 0.3;0.6e-5 4e-4; 0.06 0.2; 0.7 1.75; 0.58 0.88];
load('Range8.mat')
opts = odeset('NonNegative',1:9);
[T,~]=gsua_dpmat(odes,vars,[0 62],'Model8E','range',Range,'output',1,'opt',opts);

load('DataBello_full.mat');
ydata=DataBello.cases(519:581)';
xdata=linspace(0,length(ydata)-1,length(ydata));
ydata2=ydata;
for i=1:length(ydata)
ydata2(i)=sum(ydata(1:i));
end


my = parcluster('local');
delete(my.Jobs);
solver='lsqc';
p = parpool(str2double(getenv('SLURM_NTASKS')));
% ps = parallel.Settings;
% ps.Pool.PreferredNumWorkers=16;
% M=gsua_dmatrix(T,5000);
% [T8,J8,Y8]=gsua_sa(M,T,'SensMethod','Distance');
% save('Results8_sa.mat','T8','J8','Y8','xdata')
opt=optimoptions('lsqcurvefit','UseParallel',1,'Display','iter','MaxFunctionEvaluations',4000);
[T8,res]=gsua_pe(T,xdata,ydata2,'solver',solver,'N',100,'opt',opt);
save('Results8.mat','T8','res','xdata','ydata2')



