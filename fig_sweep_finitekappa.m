clear all
close all

addpath(genpath('utils\'))
analytical_expressions

Xsource=0.5;
R=1;

ron=0.5;
roff=.1;

kon=100;
koff=1.5;
D=5;
Tfinal=30;



nSims=5e3;
nMaxParticles=1e4;



%%
kappa=1e-3;

Nparticles1 = zeros(1,nSims);
positions = zeros(nMaxParticles,nSims);


parfor i = 1:nSims
    particles_at_end = montecarlo_telegraph_robin(Tfinal,R, Xsource,kon, koff, D, ron,roff,kappa);
    Nparticles1(i) = length(particles_at_end);
    padded_particles = [particles_at_end; NaN(nMaxParticles-length(particles_at_end),1 ) ];
    positions(:,i) = padded_particles;
end 

nBins=71;
bins = linspace(-R,R,nBins);
binned_mid = bins(2:end)*0.5+bins(1:end-1)*0.5;
binned_x = zeros(nBins-1,nSims);
bin_dx = bins(2)-bins(1);

for i = 1:nSims
binned_x(:,i) = histcounts(positions(:,i),bins)/bin_dx;
end

b1=bar(binned_mid, mean(binned_x,2),1,'EdgeAlpha',0); 
hold on;
u_analy = rescaled_uu(binned_mid,R, Xsource, D, kon, koff, kappa,ron, roff); 
binned_x1 = binned_x;

plot(binned_mid, u_analy); hold on;

%%
kappa=1;
Nparticles2 = zeros(1,nSims);
positions = zeros(nMaxParticles,nSims);


parfor i = 1:nSims
    particles_at_end = montecarlo_telegraph_robin(Tfinal,R, Xsource,kon, koff, D, ron,roff,kappa);
    Nparticles2(i) = length(particles_at_end);
    padded_particles = [particles_at_end; NaN(nMaxParticles-length(particles_at_end),1 ) ];
    positions(:,i) = padded_particles;
end 

nBins=71;
bins = linspace(-R,R,nBins);
binned_mid = bins(2:end)*0.5+bins(1:end-1)*0.5;
binned_x = zeros(nBins-1,nSims);
bin_dx = bins(2)-bins(1);

for i = 1:nSims
binned_x(:,i) = histcounts(positions(:,i),bins)/bin_dx;
end

b2=bar(binned_mid, mean(binned_x,2),1,'EdgeAlpha',0); hold on;
binned_x2 = binned_x;

u_analy = rescaled_uu(binned_mid,R, Xsource, D, kon, koff, kappa,ron, roff); 
plot(binned_mid, u_analy); hold on;



%%
kappa=1000;
Nparticles3 = zeros(1,nSims);
positions = zeros(nMaxParticles,nSims);


parfor i = 1:nSims
    particles_at_end = montecarlo_telegraph_robin(Tfinal,R, Xsource,kon, koff, D, ron,roff,kappa);
    Nparticles3(i) = length(particles_at_end);
    padded_particles = [particles_at_end; NaN(nMaxParticles-length(particles_at_end),1 ) ];
    positions(:,i) = padded_particles;
end 

nBins=71;
bins = linspace(-R,R,nBins);
binned_mid = bins(2:end)*0.5+bins(1:end-1)*0.5;
binned_x = zeros(nBins-1,nSims);
bin_dx = bins(2)-bins(1);

for i = 1:nSims
binned_x(:,i) = histcounts(positions(:,i),bins)/bin_dx;
end

b3=bar(binned_mid, mean(binned_x,2),1,'EdgeAlpha',0); hold on;
binned_x3 = binned_x;
u_analy = rescaled_uu(binned_mid,R, Xsource, D, kon, koff, kappa,ron, roff); 
plot(binned_mid, u_analy); hold on;



xlim([-R R]);
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
 set(gca,'TickLength',[0.015 0.015])


%%%%
figure;

xmax = max([Nparticles1,Nparticles2,Nparticles3]);
edges = -0.5:1:(xmax + 0.5);  % Define bin edges for integers

histogram(Nparticles1,edges,'normalization','pdf');
hold on;
histogram(Nparticles2,edges,'normalization','pdf');
histogram(Nparticles3,edges,'normalization','pdf');


x_vals = 0:1:xmax;

kappa = 1000;
ronpredict = ron_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
roffpredict = roff_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
ksynpredict  =ksyn_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);

predict_pdf = Poissbeta(ronpredict,roffpredict, ksynpredict,x_vals);
predict_pdf = predict_pdf/sum(predict_pdf);
plot(x_vals,predict_pdf);

kappa = 1;
ronpredict = ron_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
roffpredict = roff_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
ksynpredict  =ksyn_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);

predict_pdf = Poissbeta(ronpredict,roffpredict, ksynpredict,x_vals);
predict_pdf = predict_pdf/sum(predict_pdf);
plot(x_vals,predict_pdf);

kappa = 1e-3;

ronpredict = ron_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
roffpredict = roff_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
ksynpredict  =ksyn_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);

predict_pdf = Poissbeta(ronpredict,roffpredict, ksynpredict,x_vals);
predict_pdf = predict_pdf/sum(predict_pdf);
plot(x_vals,predict_pdf);


xlim([0 100]);
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
 set(gca,'TickLength',[0.015 0.015])

%%%
nParams = 13;
nSims=1e4;

paramsweep1 = logspace(-3, 3,nParams);
meanvals1 = zeros(1,nParams);
meanpredicts1 = zeros(1,nParams);

varvals1 = zeros(1,nParams);
varpredictansatz1 = zeros(1,nParams);

kldivs1 = zeros(1,nParams);
shanjens1 = zeros(1,nParams);

for j = 1:nParams
    kappa = paramsweep1(j);
    Nparticles = zeros(1,nSims);

    parfor i = 1:nSims
        particles_at_end = montecarlo_telegraph_robin(Tfinal,R, Xsource,kon, koff, D, ron,roff,kappa);
        Nparticles(i) = length(particles_at_end);
    end

    meanvals1(j) = mean(Nparticles);
    meanpredicts1(j) = mean_predict(R, Xsource, D, kon, koff, kappa, ron, roff);
    varvals1(j) = var(Nparticles);
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
shanjens1(j) = kldiv(x_vals,emp_pdf+eps,predict_pdf'+eps,'sym');  

end


figure('Position',[1000 100 800 300]);
subplot(1,3,1)
plot(paramsweep1,meanpredicts1);
hold on;
scatter(paramsweep1,meanvals1,'s');
xlim([1e-4,1e4]);
set(gca,'xscale','log')
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
 set(gca,'TickLength',[0.015 0.015])
xlim([1e-4,1e4]);
ylim([25 60])
xticks([1e-4, 1e-2, 1e0 1e2, 1e4])

subplot(1,3,2)
plot(paramsweep1,varpredictansatz1);
hold on;
scatter(paramsweep1,varvals1,'s');
xlim([1e-4,1e4]);
ylim([100 600])
set(gca,'xscale','log')
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
 set(gca,'TickLength',[0.015 0.015])
xticks([1e-4, 1e-2, 1e0 1e2, 1e4])



subplot(1,3,3)
plot(paramsweep1,varpredictansatz1./meanpredicts1);
hold on;
scatter(paramsweep1,varvals1./meanvals1,'s');
xlim([1e-4,1e4]);
%ylim([100 600])
set(gca,'xscale','log')
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
 set(gca,'TickLength',[0.015 0.015])
xticks([1e-4, 1e-2, 1e0 1e2, 1e4])





save('sweep_kappavals3')