clear all;
close all;

addpath(genpath('utils\'))
rng(1);

nData = 500;
meanR = 1;
sigma2R = 0.25;
zsigma2=0.75;

gamk = meanR^2/sigma2R;
gamthet = sigma2R/meanR;

Rvals = gamrnd(gamk,gamthet,nData,1);
zvals = betarnd(1/zsigma2,1/zsigma2,nData,1).*2.*Rvals-Rvals;

Tfinal=100;
D=0.5;

kon=250;
koff=1.5;
ron = 1; 
roff=0.25;
kappa=5;

xvals= cell(nData,1);
nVals = zeros(nData,1);


parfor i = 1:nData
    R = Rvals(i);
    Xsource = zvals(i);
    particles_at_end = montecarlo_telegraph_robin(Tfinal,R, Xsource,kon, koff, D, ron,roff,kappa);
    nVals(i) = length(particles_at_end);
    xvals{i} = particles_at_end;
end 

data.xvals = xvals;
data.Rvals = Rvals;
data.zvals = zvals;

%
%D = params(1);
%kon = params(2);
%ron=params(3);
%roff=params(4);
%kappa=params(5);
%koff=params(6);

%params = [D, kon, ron,roff, kappa, koff];

trueParams = log([kon, ron,roff, kappa, koff]);
%outlike = likelihood_robin(params, data);

%negLogL = @(logTheta) -likelihood_robin([ exp(logTheta(1)),exp(logTheta(2)),exp(logTheta(3)) ...
%    exp(logTheta(4)),  exp(logTheta(4)),  exp(logTheta(4))], data);

negLogL = @(logTheta) -likelihood_robin([ D,exp(logTheta(1)),exp(logTheta(2)) ...
    exp(logTheta(3)),  exp(logTheta(4)),  exp(logTheta(5))], data);


%theta0 = log([1 100 1 1 1 ]);
theta0 = trueParams+randn(size(trueParams));%log([100 1 1 1 1]);
options = optimset('Display','iter','TolFun',1e-8, 'TolX',1e-6,'MaxFunEvals',1e5);

logThetaEst = fminsearch(negLogL,theta0,options);
theta_est = exp(logThetaEst)


pdefactor =  @(R, Xsource) ((-1).*(D.*koff).^(1/2)+kappa+exp(1).^(2.*(D.^(-1).*koff).^(1/2).* ...
        R).*((D.*koff).^(1/2)+kappa)).^(-1).*((-1).*(D.*koff).^(1/2)+ ...
        kappa+(-1).*exp(1).^((D.^(-1).*koff).^(1/2).*(R+(-1).*Xsource)).* ...
        kappa+(-1).*exp(1).^((D.^(-1).*koff).^(1/2).*(R+Xsource)).*kappa+ ...
        exp(1).^(2.*(D.^(-1).*koff).^(1/2).*R).*((D.*koff).^(1/2)+kappa));

pdefact = pdefactor(mean(Rvals),mean(zvals));

    ron_eff = pdefact*ron/koff;
    roff_eff = pdefact*roff/koff;
    ksyn_eff = pdefact*kon/koff;
    
mean_predict = ksyn_eff*(ron_eff/(ron_eff+roff_eff));
var_predict = mean_predict + ((ron_eff*roff_eff)/((ron_eff+roff_eff)^2))*...
(ksyn_eff^2/(ron_eff+roff_eff+1));


%xmax = mean_predict*2;
%x_vals = 0:1:xmax;

%predict = Poissbeta(ron_eff,roff_eff, ksyn_eff,x_vals);
%histogram(nVals,'Normalization','pdf')

%hold on;
%plot(x_vals, predict)

%trueParams 
[profile_vals1,logthetaj_vals1]  = profile_like(negLogL, logThetaEst,1);
[profile_vals2,logthetaj_vals2]  = profile_like(negLogL, logThetaEst,2);
[profile_vals3,logthetaj_vals3]  = profile_like(negLogL, logThetaEst,3);
[profile_vals4,logthetaj_vals4]  = profile_like(negLogL, logThetaEst,4);
[profile_vals5,logthetaj_vals5]  = profile_like(negLogL, logThetaEst,5);


%plot(logthetaj_vals1,profile_vals1)
%hold on;
%plot(logthetaj_vals2,profile_vals2)
%plot(logthetaj_vals3,profile_vals3)
%plot(logthetaj_vals4,profile_vals4)
%plot(logthetaj_vals5,profile_vals5)


% 
% hold on;
% 
% D = D;
% kon = theta_est(1);
% ron=theta_est(2);
% roff=roff;
% kappa=theta_est(3);%theta_est(5);
% koff=theta_est(4);%theta_est(5);
% 
% 
% pdefactor =  @(R, Xsource) ((-1).*(D.*koff).^(1/2)+kappa+exp(1).^(2.*(D.^(-1).*koff).^(1/2).* ...
%         R).*((D.*koff).^(1/2)+kappa)).^(-1).*((-1).*(D.*koff).^(1/2)+ ...
%         kappa+(-1).*exp(1).^((D.^(-1).*koff).^(1/2).*(R+(-1).*Xsource)).* ...
%         kappa+(-1).*exp(1).^((D.^(-1).*koff).^(1/2).*(R+Xsource)).*kappa+ ...
%         exp(1).^(2.*(D.^(-1).*koff).^(1/2).*R).*((D.*koff).^(1/2)+kappa));
% 
% pdefact = pdefactor(mean(Rvals),mean(zvals));
% 
% 
%     ron_eff = pdefact*ron/koff;
%     roff_eff = pdefact*roff/koff;
%     ksyn_eff = pdefact*kon/koff;
% 
%     mean_predictnew = ksyn_eff*(ron_eff/(ron_eff+roff_eff));
% var_predictnew = mean_predictnew + ((ron_eff*roff_eff)/((ron_eff+roff_eff)^2))*...
% (ksyn_eff^2/(ron_eff+roff_eff+1));
% 
%     predict2 = Poissbeta(ron_eff,roff_eff, ksyn_eff,x_vals);
% plot(x_vals, predict2)
%figure; 
%ss = logspace(-1,1, length(profile_vals1));
%plot(2*profile_vals1); hold on
%plot(2*profile_vals2); hold on
%plot(2*profile_vals3); hold on
%plot(2*profile_vals4); hold on
%plot(2*profile_vals5); hold on

%yline(cutoff)
%ylim([0, 10*cutoff])

save('profile_run33')


cutoff= chi2inv(0.95,4)/2;
figure('Position', [1000 818 950 425]);

subplot(2,3,1);
plot(exp(logthetaj_vals1-trueParams(1)),-profile_vals1)
%set(gca,'XScale','log')
yline(-cutoff)
ylim([-20, 0])
xline(1);
pbaspect([4 3 1])
set(gca,'FontSize',11)
set(gca,'LineWidth',1.25)
box off;
%set(gca,'Xscale','log');
 set(gca,'TickLength',[0.015 0.015])
%xlim([0.5 2])


subplot(2,3,2);
plot(exp(logthetaj_vals2-trueParams(2)),-profile_vals2)
%set(gca,'XScale','log')
yline(-cutoff)
ylim([-20, 0])
xline(1);
pbaspect([4 3 1])
set(gca,'FontSize',11)
set(gca,'LineWidth',1.25)
box off;
%set(gca,'Xscale','log');
 set(gca,'TickLength',[0.015 0.015])
%xlim([0.5 2])


subplot(2,3,3);
plot(exp(logthetaj_vals3-trueParams(3)),-profile_vals3)
%set(gca,'XScale','log')
yline(-cutoff)
ylim([-20, 0])
xline(1);
pbaspect([4 3 1])
set(gca,'FontSize',11)
set(gca,'LineWidth',1.25)
box off;
%set(gca,'Xscale','log');
 set(gca,'TickLength',[0.015 0.015])
%xlim([0.5 2])


subplot(2,3,4);
plot(exp(logthetaj_vals4-trueParams(4)),-profile_vals4)
%set(gca,'XScale','log')
yline(-cutoff)
ylim([-20, 0])
xline(1);
pbaspect([4 3 1])
set(gca,'FontSize',11)
set(gca,'LineWidth',1.25)
box off;
%set(gca,'Xscale','log');
 set(gca,'TickLength',[0.015 0.015])
%xlim([0.5 2])


subplot(2,3,5);
plot(exp(logthetaj_vals5-trueParams(5)),-profile_vals5)
%set(gca,'XScale','log')
yline(-cutoff)
ylim([-20, 0])
xline(1);
pbaspect([4 3 1])
set(gca,'FontSize',11)
set(gca,'LineWidth',1.25)
box off;
%set(gca,'Xscale','log');
%xlim([0.5 2])
set(gca,'TickLength',[0.015 0.015])


% 
% subplot(2,3,2);
% plot(exp(logthetaj_vals2),profile_vals2)
% %set(gca,'XScale','log')
% xline(exp(trueParams(2)))
% 
% yline(cutoff)
% ylim([0, 20])
% 
% subplot(2,3,3);
% plot(exp(logthetaj_vals3),profile_vals3)
% %set(gca,'XScale','log')
% xline(exp(trueParams(3)))
% 
% yline(cutoff)
% ylim([0, 20])
% 
% subplot(2,3,4);
% plot(exp(logthetaj_vals4),profile_vals4)
% %set(gca,'XScale','log')
% xline(exp(trueParams(4)))
% 
% 
% yline(cutoff)
% ylim([0, 20])
% 
% 
% subplot(2,3,5);
% plot(exp(logthetaj_vals5),profile_vals5)
% %set(gca,'XScale','log')
% xline(exp(trueParams(5)))
% 
% yline(cutoff)
% ylim([0, 20])

