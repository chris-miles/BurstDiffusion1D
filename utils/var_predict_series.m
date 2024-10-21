function sumout =  var_predict_series(R, Xsource, D, kon, koff, ron, roff)
M=5000;
N=5000;
sumout = 0;
for m=1:M
    for n=1:N
        sumout = sumout + 4.*((-1)+(-1).^m).*((-1)+(-1).^n).*m.^(-1).*n.^(-1).*pi .^(-2).*( ...
            D.*(m.^2+n.^2).*pi .^2+8.*koff.*R.^2).^(-1).*(D.*m.^2.*pi .^2+4.* ...
            R.^2.*(koff+roff+ron)).^(-1).*(D.*n.^2.*pi .^2+4.*R.^2.*(koff+ ...
            roff+ron)).^(-1).*(16.*D.*(m.^2+n.^2).*pi .^2.*R.^4+128.*R.^6.*( ...
            koff+roff+ron)).*sin((1/2).*m.*pi .*R.^(-1).*(R+Xsource)).*sin(( ...
            1/2).*n.*pi .*R.^(-1).*(R+Xsource));
    end
end

sumout = sumout.*kon^2.*ron.*roff./(ron+roff).^2;

end