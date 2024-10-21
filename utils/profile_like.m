function [profile_vals,logthetaj_vals]  = profile_like(logL, logtheta0,j)

nsweep=151;
logtheta_j0 = logtheta0(j);
nparams = length(logtheta0);

logthetaj_vals = linspace(logtheta_j0-1,logtheta_j0+1,nsweep);

l0 = logL(logtheta0);
profile_vals = zeros(1,nsweep);


parfor n = 1:nsweep
    paramjn =  logthetaj_vals(n);
    optFun = @(params) feval(logL,[params(1:j-1) paramjn params(j:end)]);

    theta0 = logtheta0;
    theta0(j) = [];
    options = optimset('TolFun',1e-2, 'TolX',1e-2);
    [fittheta,logLprofnj] = fminsearch(optFun,theta0,options);
    profile_vals(n) = logLprofnj;
end

profile_vals = profile_vals-l0;