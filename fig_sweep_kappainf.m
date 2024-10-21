clear all
close all

addpath(genpath('utils\'))
analytical_expressions

kappa=1e8; % purely absorbing for now

Xsource=0;
R=3;

ron=0.5;
roff=.1;

kon=100;
koff=1.5;
D=1;
Tfinal=100;


nSims=5000;
nParams = 21;

%% zsweep
paramsweep1 = linspace(0, 0.99*R,nParams);
meanvals1 = zeros(1,nParams);
meanpredicts1 = zeros(1,nParams);

varvals1 = zeros(1,nParams);
varpredictseries1 = zeros(1,nParams);
varpredictansatz1 = zeros(1,nParams);

kldivs1 = zeros(1,nParams);
shanjens1 = zeros(1,nParams);

for j = 1:nParams
    Xsource = paramsweep1(j);

    Nparticles = zeros(1,nSims);

    parfor i = 1:nSims
        particles_at_end = montecarlo_telegraph_robin(Tfinal,R, Xsource,kon, koff, D, ron,roff,kappa);
        Nparticles(i) = length(particles_at_end);
    end

    meanvals1(j) = mean(Nparticles);
    meanpredicts1(j) = mean_predict(R, Xsource, D, kon, koff, kappa, ron, roff);
    varvals1(j) = var(Nparticles);
    varpredictseries1(j) = var_predict_series(R, Xsource, D, kon, koff,  ron, roff) + mean_predict(R, Xsource, D, kon, koff, kappa, ron, roff);
    varpredictansatz1(j) = var_predict(R, Xsource, D, kon, koff, kappa,  ron, roff);


    xmax = max(Nparticles);
    x_vals = 0:1:xmax;

    ronpredict = ron_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
    roffpredict = roff_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
    ksynpredict  =ksyn_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);

    predict_pdf = Poissbeta(ronpredict,roffpredict, ksynpredict,x_vals);
    predict_pdf = predict_pdf/sum(predict_pdf);
    edges = -0.5:1:(xmax + 0.5);  % Define bin edges for integers

    counts = histcounts(Nparticles, edges);  % Get counts for each bin
    emp_pdf = counts / length(Nparticles);

    kldivs1(j) = kldiv(x_vals,emp_pdf+eps,predict_pdf'+eps);
    shanjens1(j) = kldiv(x_vals,emp_pdf+eps,predict_pdf'+eps,'js');

end



%% gamma sweep
Xsource=0;
paramsweep2 = logspace(-1, 1, nParams);
meanvals2 = zeros(1,nParams);
meanpredicts2 = zeros(1,nParams);

varvals2 = zeros(1,nParams);
varpredictseries2 = zeros(1,nParams);
varpredictansatz2 = zeros(1,nParams);

kldivs2 = zeros(1,nParams);
shanjens2 = zeros(1,nParams);

for j = 1:nParams
    koff = paramsweep2(j);
    Nparticles = zeros(1,nSims);

    parfor i = 1:nSims
        particles_at_end = montecarlo_telegraph_robin(Tfinal,R, Xsource,kon, koff, D, ron,roff,kappa);
        Nparticles(i) = length(particles_at_end);
    end

    meanvals2(j) = mean(Nparticles);
    meanpredicts2(j) = mean_predict(R, Xsource, D, kon, koff, kappa, ron, roff);
    varvals2(j) = var(Nparticles);
    varpredictseries2(j) = var_predict_series(R, Xsource, D, kon, koff,  ron, roff) + mean_predict(R, Xsource, D, kon, koff, kappa, ron, roff);
    varpredictansatz2(j) = var_predict(R, Xsource, D, kon, koff, kappa,  ron, roff);

    xmax = max(Nparticles);
    x_vals = 0:1:xmax;

    ronpredict = ron_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
    roffpredict = roff_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
    ksynpredict  =ksyn_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);

    predict_pdf = Poissbeta(ronpredict,roffpredict, ksynpredict,x_vals);
    predict_pdf = predict_pdf/sum(predict_pdf);
    edges = -0.5:1:(xmax + 0.5);  % Define bin edges for integers

    counts = histcounts(Nparticles, edges);  % Get counts for each bin
    emp_pdf = counts / length(Nparticles);

    kldivs2(j) = kldiv(x_vals,emp_pdf+eps,predict_pdf'+eps);
    shanjens2(j) = kldiv(x_vals,emp_pdf+eps,predict_pdf'+eps,'js');

end




%% D sweep
koff=1.5;
paramsweep3 = logspace(-1, 2, nParams);
meanvals3 = zeros(1,nParams);
meanpredicts3 = zeros(1,nParams);

varvals3 = zeros(1,nParams);
varpredictseries3 = zeros(1,nParams);
varpredictansatz3 = zeros(1,nParams);


kldivs3 = zeros(1,nParams);
shanjens3 = zeros(1,nParams);

for j = 1:nParams
    D = paramsweep3(j);

    Nparticles = zeros(1,nSims);

    parfor i = 1:nSims
        particles_at_end = montecarlo_telegraph_robin(Tfinal,R, Xsource,kon, koff, D, ron,roff,kappa);
        Nparticles(i) = length(particles_at_end);
    end

    meanvals3(j) = mean(Nparticles);
    meanpredicts3(j) = mean_predict(R, Xsource, D, kon, koff, kappa, ron, roff);
    varvals3(j) = var(Nparticles);
    varpredictseries3(j) = var_predict_series(R, Xsource, D, kon, koff,  ron, roff) + mean_predict(R, Xsource, D, kon, koff, kappa, ron, roff);
    varpredictansatz3(j) = var_predict(R, Xsource, D, kon, koff, kappa,  ron, roff);

    xmax = max(Nparticles);
    x_vals = 0:1:xmax;

    ronpredict = ron_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
    roffpredict = roff_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
    ksynpredict  =ksyn_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);

    predict_pdf = Poissbeta(ronpredict,roffpredict, ksynpredict,x_vals);
    predict_pdf = predict_pdf/sum(predict_pdf);
    edges = -0.5:1:(xmax + 0.5);  % Define bin edges for integers

    counts = histcounts(Nparticles, edges);  % Get counts for each bin
    emp_pdf = counts / length(Nparticles);

    kldivs3(j) = kldiv(x_vals,emp_pdf+eps,predict_pdf'+eps);
    shanjens3(j) = kldiv(x_vals,emp_pdf+eps,predict_pdf'+eps,'js');


end


%% ron sweep
D=1;
paramsweep4 = logspace(-1, 1, nParams);

meanvals4 = zeros(1,nParams);
meanpredicts4 = zeros(1,nParams);

varvals4 = zeros(1,nParams);
varpredictseries4 = zeros(1,nParams);
varpredictansatz4 = zeros(1,nParams);


kldivs4 = zeros(1,nParams);
shanjens4 = zeros(1,nParams);


for j = 1:nParams
    ron = paramsweep4(j);
    Nparticles = zeros(1,nSims);

    parfor i = 1:nSims
        particles_at_end = montecarlo_telegraph_robin(Tfinal,R, Xsource,kon, koff, D, ron,roff,kappa);
        Nparticles(i) = length(particles_at_end);
    end

    meanvals4(j) = mean(Nparticles);
    meanpredicts4(j) = mean_predict(R, Xsource, D, kon, koff, kappa, ron, roff);
    varvals4(j) = var(Nparticles);
    varpredictseries4(j) = var_predict_series(R, Xsource, D, kon, koff,  ron, roff) + mean_predict(R, Xsource, D, kon, koff, kappa, ron, roff);
    varpredictansatz4(j) = var_predict(R, Xsource, D, kon, koff, kappa,  ron, roff);

    xmax = max(Nparticles);
    x_vals = 0:1:xmax;

    ronpredict = ron_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
    roffpredict = roff_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
    ksynpredict  =ksyn_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);

    predict_pdf = Poissbeta(ronpredict,roffpredict, ksynpredict,x_vals);
    predict_pdf = predict_pdf/sum(predict_pdf);
    edges = -0.5:1:(xmax + 0.5);  % Define bin edges for integers

    counts = histcounts(Nparticles, edges);  % Get counts for each bin
    emp_pdf = counts / length(Nparticles);

    kldivs4(j) = kldiv(x_vals,emp_pdf+eps,predict_pdf'+eps);
    shanjens4(j) = kldiv(x_vals,emp_pdf+eps,predict_pdf'+eps,'js');


end


figure('position',[1 1 850 650]);
subplot(2,2,1);
plot(paramsweep1, meanpredicts1)
hold on;
scatter(paramsweep1, meanvals1,'s');
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
set(gca,'TickLength',[0.015 0.015])

subplot(2,2,2);
plot(paramsweep2, meanpredicts2)
hold on;
scatter(paramsweep2, meanvals2,'s');
set(gca,'XScale','log')

pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
set(gca,'TickLength',[0.015 0.015])


subplot(2,2,3);
plot(paramsweep3, meanpredicts3)
hold on;
scatter(paramsweep3, meanvals3,'s');
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
set(gca,'TickLength',[0.015 0.015])
set(gca,'XScale','log')


subplot(2,2,4);
plot(paramsweep4, meanpredicts4)
hold on;
scatter(paramsweep4, meanvals4,'s');
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
set(gca,'XScale','log')
set(gca,'TickLength',[0.015 0.015])

figure('position',[1 1 850 650]);
subplot(2,2,1);
plot(paramsweep1, varpredictseries1)
hold on;
plot(paramsweep1, varpredictansatz1,'--')
scatter(paramsweep1, varvals1,'s');
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
set(gca,'TickLength',[0.015 0.015])

subplot(2,2,2);
plot(paramsweep2, varpredictseries2)
hold on;
plot(paramsweep2, varpredictansatz2,'--')
scatter(paramsweep2, varvals2,'s');
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
set(gca,'TickLength',[0.015 0.015])
set(gca,'XScale','log')

subplot(2,2,3);
plot(paramsweep3, varpredictseries3)
hold on;
plot(paramsweep3, varpredictansatz3,'--')
scatter(paramsweep3, varvals3,'s');
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
set(gca,'TickLength',[0.015 0.015])
set(gca,'XScale','log')

subplot(2,2,4);
plot(paramsweep4, varpredictseries4)
hold on;
plot(paramsweep4, varpredictansatz4,'--')
scatter(paramsweep4, varvals4,'s');
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
set(gca,'TickLength',[0.015 0.015])
set(gca,'XScale','log')


figure('position',[1 1 850 650]);
subplot(2,2,1);
%plot(paramsweep1, kldivs1)
%hold on;
plot(paramsweep1, shanjens1)
pbaspect([4 3 1])
ylim([0 1]);
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
set(gca,'TickLength',[0.015 0.015])

subplot(2,2,2);
%plot(paramsweep1, kldivs2)
%hold on;
plot(paramsweep2, shanjens2)
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
set(gca,'TickLength',[0.015 0.015])
ylim([0 1]);
set(gca,'XScale','log')


subplot(2,2,3);
%plot(paramsweep1, kldivs3)
%hold on;
plot(paramsweep3, shanjens3)
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
set(gca,'TickLength',[0.015 0.015])
ylim([0 1]);
set(gca,'XScale','log')


subplot(2,2,4);
%plot(paramsweep1, kldivs4)
%hold on;
plot(paramsweep4, shanjens4)
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
set(gca,'TickLength',[0.015 0.015])
ylim([0 1]);
set(gca,'XScale','log')




save('paramsweeps_N_kappa_inf2')