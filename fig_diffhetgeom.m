clear all
close all

addpath(genpath('utils\'))


source_loc = [0,-0.5];
Dfast=.1;
Dslow = Dfast/100;

%Dfunc = @(x,y) Dslow + (Dfast-Dslow).*1./(1+exp(-dsteep*(sqrt(x.^2+y.^2)-rD)));
Dfunc = @(x,y) Dslow+(Dfast-Dslow)*(cos(3*pi*x) + cos(3*pi*y)  + 2)/4;

kon=100;
koff=1;
ron = 0.75; 
roff=0.25;
R=1;
Xsource=1;

Tfinal=30;

nSims=5000;

%https://stackoverflow.com/questions/54828017/c-create-random-shaped-blob-objects
rng(1);

npts=1000;
thetavals = linspace(0,2*pi,npts);

ss=5;
N=8;
amps = ss*rand(1,N)/(2*N);
phases = rand(1,N)*2*pi;

radius=1;
for j =0:N-1
    radius = radius+ amps(j+1).* cos((j+1)*thetavals + phases(j+1));
end




bdry = [cos(thetavals).*radius;sin(thetavals).*radius];

spacing=0.01;
[X, Y] = meshgrid(min(bdry(1,:)):spacing:max(bdry(1,:)),min(bdry(2,:)):spacing:max(bdry(2,:)));

mask = zeros(size(X));
Dmask = zeros(size(X));

for i=1:size(X,1)
    for j = 1:size(Y,2)
        pt=[X(i,j),Y(i,j)];
        mask(i,j) = insidepoly(pt(1),pt(2),bdry(1,:),bdry(2,:));
        if mask(i,j)==0
            Dmask(i,j) = NaN;
        else
            Dmask(i,j) = Dfunc(X(i,j),Y(i,j));
        end
    end
end

Dmaskcopy = Dmask;

Dmaskcopy(isnan(Dmaskcopy))=0;
[gradDx,gradDy] = gradient(Dmaskcopy,spacing,spacing);
gradDx_interp = griddedInterpolant(X',Y', gradDx');
gradDy_interp = griddedInterpolant(X',Y', gradDy');
gradD = @(x,y) [gradDx_interp(x,y), gradDy_interp(x,y)];

mask_flat = find(mask);
mask_x = X(mask_flat);
mask_y = Y(mask_flat);
pts_inside = [mask_x,mask_y];

shp = alphaShape(pts_inside,'HoleThreshold',1);
p = shp.Points';

tri = alphaTriangulation(shp)';
model = createpde(1);

gg=geometryFromMesh(model,p,tri);

mesh=generateMesh(model,'Hmax',.05);%,VertexID2,7},'Hmax',15)


sigma=.1;
ff = @(location,state) (3/(pi*sigma^2))*(1 - 1.0.*sqrt(sigma^(-2)*(location.x-source_loc(1)).^2+sigma^(-2)*(location.y-source_loc(2)).^2)).* ...
    1.0.*(sqrt((location.x-source_loc(1)).^2+(location.y-source_loc(2)).^2)<(sigma));

f = @(location,state) ff(location,state);

nEdges = model.Geometry.NumEdges;
applyBoundaryCondition(model,"dirichlet","Edge",(1:nEdges),"u",0);

m=0;
d=0;
a=koff;

ctot = @(location,state) [Dfunc(location.x,location.y) + (state.u).*gradDx_interp(location.x,location.y).*(state.ux+1e-24).^(-1);
    Dfunc(location.x,location.y) + (state.u).*gradDy_interp(location.x,location.y).*(state.uy+1e-24).^(-1)];

specifyCoefficients(model,'m',m,'d',d,'c',ctot,'a',a,'f',f);

model.SolverOptions.ReportStatistics = 'on';
%model.SolverOptions.ResidualTolerance = 1e-3;
%model.SolverOptions.MinStep = 1e-32;
%model.SolverOptions.ResidualNorm = 1;
%model.SolverOptions.MaxIterations = 100;

results = solvepde(model);
nmean = compute_volume(results);


sol = results.NodalSolution;
pdeplot(model,'XYData',sol);
axis equal;

set(gca,'ColorScale','log');
colormap turbo
clim([1e-3 10])


Nparticles = zeros(1,nSims);
nMaxParticles=1000;
positions = zeros(nSims,nMaxParticles,2);
positions(:)=NaN;

parfor i = 1:nSims
    particles_at_end = montecarlo_telegraph_2d_diffhet(Tfinal,bdry, source_loc, kon, koff, Dfunc, ron,roff);
    Nparticles(i) = length(particles_at_end);
    padded_particles = [particles_at_end; NaN(nMaxParticles-size(particles_at_end,1),2) ];
    positions(i,:,:) = padded_particles;
end

flattened_pos = reshape(positions,nSims*nMaxParticles,2);
flattened_pos = rmmissing(flattened_pos,1);



figure;

pdefactor = koff*nmean;

mean_pde = (ron/(ron+roff))*(kon/koff)*pdefactor;  %(-kon*(sech(R*sqrt(koff/D))-1)/koff);
mean_nopde =  (ron/(ron+roff))*(kon/koff);

scale_factor = mean_pde/mean_nopde;

xmax = mean_pde*3;
x_vals = 0:1:xmax;

ron_eff = scale_factor*ron/koff;
roff_eff = scale_factor*roff/koff;
ksyn_eff = scale_factor*kon/koff;

mean_predict = ksyn_eff*(ron_eff/(ron_eff+roff_eff));
var_predict = mean_predict + ((ron_eff*roff_eff)/((ron_eff+roff_eff)^2))*...
    (ksyn_eff^2/(ron_eff+roff_eff+1));

predict = Poissbeta(ron_eff,roff_eff, ksyn_eff,x_vals);
%edges = -0.5:1:(xmax + 0.5);  % Define bin edges for integers
histogram(Nparticles,'Normalization','pdf')

hold on;
plot(x_vals, predict)
pbaspect([4 3 1])
set(gca,'FontSize',13)
set(gca,'LineWidth',1.25)
box off;
 set(gca,'TickLength',[0.015 0.015])
 xlim([0 150])

 figure;
 scatter_kde(flattened_pos(1:1e4,1),flattened_pos(1:1e4,2),'MarkerSize',10,'filled');
 hold on;
 plot(bdry(1,:),bdry(2,:));
 axis equal