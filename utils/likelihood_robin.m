function out_tot =  likelihood_robin(params, data)


positions = data.xvals;
nData = length(positions);
Rvals = data.Rvals;
zvals = data.zvals;


D = params(1);
kon = params(2);
ron=params(3);
roff=params(4);
kappa=params(5);
koff=params(6);

pvals = zeros(1,nData);


pdefactor =  @(R, Xsource) ((-1).*(D.*koff).^(1/2)+kappa+exp(1).^(2.*(D.^(-1).*koff).^(1/2).* ...
        R).*((D.*koff).^(1/2)+kappa)).^(-1).*((-1).*(D.*koff).^(1/2)+ ...
        kappa+(-1).*exp(1).^((D.^(-1).*koff).^(1/2).*(R+(-1).*Xsource)).* ...
        kappa+(-1).*exp(1).^((D.^(-1).*koff).^(1/2).*(R+Xsource)).*kappa+ ...
        exp(1).^(2.*(D.^(-1).*koff).^(1/2).*R).*((D.*koff).^(1/2)+kappa));


    uu = @(x, R, Xsource) (1/2).*exp(1).^((-1).*(D.^(-1).*koff).^(1/2).*(2.*R+x+Xsource)).*( ...
        D.*koff).^(-1/2).*((D.*koff).^(1/2)+kappa).^(-1).*(D.*((-1)+exp(1) ...
        .^(4.*(D.^(-1).*koff).^(1/2).*R)).*koff+kappa.*(2.*(D.*koff).^( ...
        1/2)+(-1).*kappa+exp(1).^(4.*(D.^(-1).*koff).^(1/2).*R).*(2.*(D.* ...
        koff).^(1/2)+kappa))).^(-1).*((-1).*kappa.^2.*(exp(1).^(2.*(D.^( ...
        -1).*koff).^(1/2).*(R+Xsource)).*((D.*koff).^(1/2)+(-1).*kappa)+ ...
        exp(1).^(4.*(D.^(-1).*koff).^(1/2).*R).*((D.*koff).^(1/2)+kappa)+ ...
        exp(1).^(2.*(D.^(-1).*koff).^(1/2).*(2.*R+x+Xsource)).*((D.*koff) ...
        .^(1/2)+kappa)+(-1).*exp(1).^(2.*(D.^(-1).*koff).^(1/2).*(3.*R+x)) ...
        .*(3.*(D.*koff).^(1/2)+kappa))+D.*koff.*(exp(1).^(2.*(D.^(-1).* ...
        koff).^(1/2).*(R+Xsource)).*((D.*koff).^(1/2)+(-1).*kappa)+exp(1) ...
        .^(4.*(D.^(-1).*koff).^(1/2).*R).*((D.*koff).^(1/2)+kappa)+exp(1) ...
        .^(2.*(D.^(-1).*koff).^(1/2).*(2.*R+x+Xsource)).*((D.*koff).^(1/2) ...
        +kappa)+exp(1).^(2.*(D.^(-1).*koff).^(1/2).*(3.*R+x)).*((D.*koff) ...
        .^(1/2)+3.*kappa))+(-1).*(exp(1).^(2.*(D.^(-1).*koff).^(1/2).*(R+ ...
        x))+(-1).*exp(1).^(2.*(D.^(-1).*koff).^(1/2).*(R+Xsource))).*( ...
        kappa.^2.*((D.*koff).^(1/2)+(-1).*kappa+exp(1).^(4.*(D.^(-1).* ...
        koff).^(1/2).*R).*(3.*(D.*koff).^(1/2)+kappa))+D.*koff.*((-1).*( ...
        D.*koff).^(1/2)+kappa+exp(1).^(4.*(D.^(-1).*koff).^(1/2).*R).*(( ...
        D.*koff).^(1/2)+3.*kappa))).*heaviside(x+(-1).*Xsource));

parfor n = 1:nData
    pos_n = positions{n};
    nn = length(pos_n);
    R = Rvals(n);
    Xsource = zvals(n);

    pdefact_n = pdefactor(R,Xsource);
    
    %rescaled_uu = @(x) koff*uu(x)/pdefactor(R,Xsource);

    ron_eff = pdefact_n*ron/koff;
    roff_eff = pdefact_n*roff/koff;
    ksyn_eff = pdefact_n*kon/koff;



%mean_predict = ksyn_eff*(ron_eff/(ron_eff+roff_eff));
%var_predict = mean_predict + ((ron_eff*roff_eff)/((ron_eff+roff_eff)^2))*...
%(ksyn_eff^2/(ron_eff+roff_eff+1));

    pn = log(Poissbeta(ron_eff,roff_eff, ksyn_eff,nn));
    
    if nn >0
        pxvals = koff*uu(pos_n, R,Xsource)/pdefact_n;
        px = nansum(log(pxvals)); % mean seems to work better
    else
      px =0 
    end 
   % for j = 1:nn
   %     xx = pos_n(j);
    %    px = px + log(rescaled_uu(xx));
    %end

    pvals(n) = px+pn;

end

out_tot = nanmean(pvals)*nData;

end