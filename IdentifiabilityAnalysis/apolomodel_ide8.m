addpath Tool/

load('best_outputs.mat');
ywork=y8;

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
load('Range8.mat')
opts = odeset('NonNegative',1:9);
[T8,~]=gsua_dpmat(odes,vars,[0 62],'Model8E','range',Range,'output',1,'opt',opts);

xdata=linspace(0,length(ywork)-1,length(ywork));
T8.Estlsqc=est8;

my = parcluster('local');
delete(my.Jobs);
solver='lsqc';
p = parpool(str2double(getenv('SLURM_NTASKS')));
% ps = parallel.Settings;
% ps.Pool.PreferredNumWorkers=16;
opt=optimoptions('lsqcurvefit','UseParallel',true,'MaxFunctionEvaluations',5000,'MaxIterations',400,'Display','iter');
n=1000;
exp=[1,5,10,20,40];
len=length(exp);
res=zeros(size(T8,1),n,len);
signal=zeros(n,length(ywork),len);
resT=T8(:,'Range');


for j=1:len-2
pd=makedist('normal','mu',0,'sigma',exp(j)/100);
ydata=ywork+ywork.*random(pd,n,length(ywork));
ydata2=ydata;
for i=1:size(ydata,2)
ydata2(:,i)=sum(ydata(:,1:i),2);
end
signal(:,:,j)=ydata2;

for i=1:n
    clc
    disp(strcat(num2str((i+((j-1)*n))/(n*len)),'%'))
    [T,~]=gsua_pe(T8,xdata,signal(i,:,j),'solver','lsqc','ipoint',T8.Estlsqc(:,1)','opt',opt);
    res(:,i,j)=T.Estlsqc;
end
are=100*sum(abs(res(:,:,j)-T8.Estlsqc(:,1)),2)./T8.Estlsqc(:,1)/n;
eval(strcat('resT.are',num2str(exp(j)),'=round(are,1);'));
save('Ident8to8.mat','resT','res','signal')
end
