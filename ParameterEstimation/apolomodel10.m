addpath Tool/
syms E(t) L(t) P(t) Ms(t) Me(t) Mi(t) Hs(t) He(t) Hi(t) Hr(t) Hit(t)

syms alpha delta C gamma_e mu_e gamma_l mu_l gamma_p mu_p f beta_m mu_m theta_m mu_h beta_h theta_h...
    gamma_h

%pulse = piecewise(t0_c<=t<=t0_c+Delta_tc,A_c,0);
ode1 = diff(E) == delta*(1-E/C)*(Ms+Me+Mi) - (gamma_e+mu_e)*E;
ode2 = diff(L) == gamma_e*E - (gamma_l+mu_l)*L;
ode3 = diff(P) == gamma_l*L - (gamma_p+mu_p)*P;
ode4 = diff(Ms) == f*gamma_p*P - beta_m*Hi*Ms/(Hs+He+Hi+Hr) - (mu_m)*Ms;
ode5 = diff(Me) == beta_m*Hi*Ms/(Hs+He+Hi+Hr) - (theta_m+(mu_m))*Me;
ode6 = diff(Mi) == theta_m*Me - (mu_m)*Mi;
ode7 = diff(Hs) == -beta_h*Mi*Hs/(Ms+Me+Mi) + (He+Hi+Hr)*mu_h;
ode8 = diff(He) == beta_h*Mi*Hs/(Ms+Me+Mi) - (theta_h+mu_h)*He;
ode9 = diff(Hi) == theta_h*He - (gamma_h+mu_h)*Hi;
ode10 = diff(Hr) == gamma_h*Hi - mu_h*Hr;
ode11 = diff(Hit) == theta_h*He;
odes=[ode1; ode2 ;ode3; ode4; ode5 ;ode6; ode7; ode8; ode9; ode10; ode11];
vars=[Hit Hi E L P Ms Me Mi Hs He Hr];

load('Range10.mat')
% range=[0 3e4; 0 3e4; 0 3e4; 1e3 8e6; 0 100; 0 100; 24.44e4 32.17e4; 0 100; 6 24; 8.14e4 15.88e4;64e2 34e4; 1 1.5; 0 4; 0 4; 20 200; 0.3 0.7;...
%     0.6 2.3; 0.3 2; 0.05 0.5; 0.1 1; 0 1.3; 1e-05 9e-04; 0.07 3.22; 0 0.9; 0 1; 0.5 1.8; 0.4 0.9];
opts = odeset('NonNegative',1:11);
[T,~]=gsua_dpmat(odes,vars,[0 62],'10p','output',1,'opt',opts,'range',Range);
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
% M=gsua_dmatrix(T,5000);
% [T10,J10,Y10]=gsua_sa(M,T,'SensMethod','Distance');
% ps = parallel.Settings;
% ps.Pool.PreferredNumWorkers=16;
opt=optimoptions('lsqcurvefit','UseParallel',1,'Display','iter','MaxFunctionEvaluations',4000);
[T10,res]=gsua_pe(T,xdata,ydata2,'solver',solver,'N',100,'opt',opt);
% save('Results10_sa.mat','T10','J10','Y10','xdata')
save('Results10.mat','T10','res','xdata','ydata2')

% inter=intersect(T.Properties.RowNames,T2.Properties.RowNames)
% T(inter,'Range')=T2(inter,'Range')