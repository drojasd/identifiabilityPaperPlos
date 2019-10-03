addpath Tool/

load('best_outputs.mat');
ywork=y7_fix;

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
[T7,~]=gsua_dpmat(odes,vars,[0 62],'7m2','output',1,'opt',opts,'Range',Range);


xdata=linspace(0,length(ywork)-1,length(ywork));
T7.Estlsqc=est7_fix;

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
res=zeros(size(T7,1),n,len);
signal=zeros(n,length(ywork),len);
resT=T7(:,'Range');


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
    [T,~]=gsua_pe(T7,xdata,signal(i,:,j),'solver','lsqc','ipoint',T7.Estlsqc(:,1)','opt',opt);
    res(:,i,j)=T.Estlsqc;
end
are=100*sum(abs(res(:,:,j)-T7.Estlsqc(:,1)),2)./T7.Estlsqc(:,1)/n;
eval(strcat('resT.are',num2str(exp(j)),'=round(are,1);'));
save('Ident7to7_fix.mat','resT','res','signal')
end
