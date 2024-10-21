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
Tfinal=30;



%%
nSims=1e4;
nMaxParticles=1e4;

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

b=bar(binned_mid, mean(binned_x,2),1,'EdgeAlpha',0); hold on;
u_analy = rescaled_uu(binned_mid,R, Xsource, D, kon, koff, kappa,ron, roff); 
plot(binned_mid, u_analy); hold on;

%%
Xsource=1;
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

b=bar(binned_mid, mean(binned_x,2),1,'EdgeAlpha',0); hold on;
u_analy = rescaled_uu(binned_mid,R, Xsource, D, kon, koff, kappa,ron, roff); 
plot(binned_mid, u_analy); hold on;



%%
Xsource=2;
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

b=bar(binned_mid, mean(binned_x,2),1,'EdgeAlpha',0); hold on;
u_analy = rescaled_uu(binned_mid,R, Xsource, D, kon, koff, kappa,ron, roff); 
plot(binned_mid, u_analy); hold on;


%%
Xsource=2.5;
Nparticles4 = zeros(1,nSims);
positions = zeros(nMaxParticles,nSims);


parfor i = 1:nSims
    particles_at_end = montecarlo_telegraph_robin(Tfinal,R, Xsource,kon, koff, D, ron,roff,kappa);
    Nparticles4(i) = length(particles_at_end);
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

b=bar(binned_mid, mean(binned_x,2),1,'EdgeAlpha',0); hold on;
u_analy = rescaled_uu(binned_mid,R, Xsource, D, kon, koff, kappa,ron, roff); 
plot(binned_mid, u_analy); hold on;
pbaspect([5 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
 set(gca,'TickLength',[0.015 0.015])


 figure;

xmax = max([Nparticles1,Nparticles2,Nparticles3,Nparticles4]);
edges = -0.5:1:(xmax + 0.5);  % Define bin edges for integers

histogram(Nparticles1,edges,'normalization','pdf');
hold on;
histogram(Nparticles2,edges,'normalization','pdf');
histogram(Nparticles3,edges,'normalization','pdf');
histogram(Nparticles4,edges,'normalization','pdf');


x_vals = 0:1:xmax;

Xsource = 0;
ronpredict = ron_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
roffpredict = roff_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
ksynpredict  =ksyn_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);

predict_pdf = Poissbeta(ronpredict,roffpredict, ksynpredict,x_vals);
predict_pdf = predict_pdf/sum(predict_pdf);
plot(x_vals,predict_pdf);

Xsource = 1;
ronpredict = ron_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
roffpredict = roff_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
ksynpredict  =ksyn_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);

predict_pdf = Poissbeta(ronpredict,roffpredict, ksynpredict,x_vals);
predict_pdf = predict_pdf/sum(predict_pdf);
plot(x_vals,predict_pdf);

Xsource = 2;
ronpredict = ron_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
roffpredict = roff_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
ksynpredict  =ksyn_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);

predict_pdf = Poissbeta(ronpredict,roffpredict, ksynpredict,x_vals);
predict_pdf = predict_pdf/sum(predict_pdf);
plot(x_vals,predict_pdf);

Xsource = 2.5;
ronpredict = ron_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
roffpredict = roff_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);
ksynpredict  =ksyn_eff(R, Xsource, D, kon, koff, kappa,  ron, roff);

predict_pdf = Poissbeta(ronpredict,roffpredict, ksynpredict,x_vals);
predict_pdf = predict_pdf/sum(predict_pdf);
plot(x_vals,predict_pdf);

pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
 set(gca,'TickLength',[0.015 0.015])
xlim([0 100])