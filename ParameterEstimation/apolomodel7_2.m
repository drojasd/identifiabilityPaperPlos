addpath Tool/
syms  Ms(t) Me(t) Mi(t) Hs(t) He(t) Hi(t) Hr(t) Hit(t)

syms lambda beta_m mu_m theta_m mu_h beta_h theta_h...
    gamma_h

%pulse(t)=1;
H=Hs+He+Hi+Hr;
M=Ms+Me+Mi;
ode1 = diff(Ms) == lambda - beta_m*Hi*Ms/H - (mu_m)*Ms;
ode2 = diff(Me) == beta_m*Hi*Ms/H - (theta_m+mu_m)*Me;
ode3 = diff(Mi) == theta_m*Me - mu_m*Mi;
ode4 = diff(Hs) == -beta_h*Mi*Hs/M + (He+Hi+Hr)*mu_h;
ode5 = diff(He) == beta_h*Mi*Hs/M - (theta_h+mu_h)*He;
ode6 = diff(Hi) == theta_h*He - (gamma_h+mu_h)*Hi;
ode7 = diff(Hr) == gamma_h*Hi - mu_h*Hr;
ode8 = diff(Hit) == theta_h*He;
odes=[ode1; ode2 ;ode3; ode4; ode5 ;ode6; ode7; ode8];
vars=[Hit Hi Me Hr Hs He Ms Mi];
opts = odeset('NonNegative',1:8);
load Range7_2.mat
[T,~]=gsua_dpmat(odes,vars,[0 62],'7m2','output',1,'opt',opts,'Range',Range);
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
% [T7,J7,Y7]=gsua_sa(M,T,'SensMethod','Distance');
% save('Results7_sa.mat','T7','J7','Y7','xdata')
opt=optimoptions('lsqcurvefit','UseParallel',1,'Display','iter','MaxFunctionEvaluations',4000);
[T7,res]=gsua_pe(T,xdata,ydata2,'solver',solver,'N',1000,'opt',opt);
save('Results7_ref.mat','T7','res','xdata','ydata2')