addpath Tool/

load('best_outputs.mat');
ywork=y10;

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
opts = odeset('NonNegative',1:9);
[T10,~]=gsua_dpmat(odes,vars,[0 62],'Model8E','range',Range,'output',1,'opt',opts);

xdata=linspace(0,length(ywork)-1,length(ywork));
T10.Estlsqc=est10;

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
res=zeros(size(T10,1),n,len);
signal=zeros(n,length(ywork),len);
resT=T10(:,'Range');


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
    [T,~]=gsua_pe(T10,xdata,signal(i,:,j),'solver','lsqc','ipoint',T10.Estlsqc(:,1)','opt',opt);
    res(:,i,j)=T.Estlsqc;
end
are=100*sum(abs(res(:,:,j)-T10.Estlsqc(:,1)),2)./T10.Estlsqc(:,1)/n;
eval(strcat('resT.are',num2str(exp(j)),'=round(are,1);'));
save('Ident10to10.mat','resT','res','signal')
end
