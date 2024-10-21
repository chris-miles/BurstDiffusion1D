clear all
close all

addpath(genpath('utils'))
analytical_expressions

kon=100;
ron=0.5;
roff=0.1;
D=1;
R=1;
koff=0.75;
kappa=10000;


figure('Position',[100, 100, 1700, 350])
subplot(1,4,1);
zvals = linspace(-1,1,25);

vv = var_predict(R, zvals, D, kon, koff, kappa, ron, roff);
mm = mean_predict(R, zvals, D, kon, koff, kappa, ron, roff);

fano_spat = vv./mm;
fano_nonspat = 1+(roff*kon)./((ron+roff)*(ron+roff+koff));
plot(zvals, fano_spat/fano_nonspat); hold on;
ylim([0 1]);

kappa=1;
vv = var_predict(R, zvals, D, kon, koff, kappa, ron, roff);
mm = mean_predict(R, zvals, D, kon, koff, kappa, ron, roff);

fano_spat = vv./mm;
plot(zvals, fano_spat/fano_nonspat); hold on;
ylim([0 1]);

 kappa=1e-3;
 vv = var_predict(R, zvals, D, kon, koff, kappa, ron, roff);
 mm = mean_predict(R, zvals, D, kon, koff, kappa, ron, roff);
 
 fano_spat = vv./mm;
 plot(zvals, fano_spat/fano_nonspat); hold on;
 ylim([0 1]);

ylim([0 1]);
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
 set(gca,'TickLength',[0.015 0.015])
 %%
z=0;
kappa=1000;

subplot(1,4,2)
koffvals = logspace(-1,2,25);

vv = var_predict(R, z, D, kon, koffvals, kappa, ron, roff);
mm = mean_predict(R, z, D, kon, koffvals, kappa, ron, roff);

fano_spat = vv./mm;
fano_nonspat = 1+(roff*kon)./((ron+roff).*(ron+roff+koffvals));
plot(koffvals, fano_spat./fano_nonspat); hold on;
ylim([0 1]);
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
set(gca,'Xscale','log');
 set(gca,'TickLength',[0.015 0.015])
%%
z=0;
koff=1.5;

subplot(1,4,3)

omegavals = logspace(-1,2,25);
p=0.75;
ron=p*omegavals;
roff=omegavals-p*omegavals;

vv = var_predict(R, z, D, kon, koff, kappa, ron, roff);
mm = mean_predict(R, z, D, kon, koff, kappa, ron, roff);

fano_spat = vv./mm;
fano_nonspat = 1+(roff*kon)./((ron+roff).*(ron+roff+koff));
plot(omegavals, fano_spat./fano_nonspat); hold on;
ylim([0 1]);
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
set(gca,'Xscale','log');
 set(gca,'TickLength',[0.015 0.015])
%%

subplot(1,4,4)

pvals = linspace(0,1,251);
omega=0.1;
ron=pvals*omega;
roff=omega-pvals*omega;

vv = var_predict(R, z, D, kon, koff, kappa, ron, roff);
mm = mean_predict(R, z, D, kon, koff, kappa, ron, roff);

fano_spat = vv./mm;
fano_nonspat = 1+(roff*kon)./((ron+roff).*(ron+roff+koff));
plot(pvals, fano_spat./fano_nonspat); hold on;
ylim([0 1]);
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
 set(gca,'TickLength',[0.015 0.015])
