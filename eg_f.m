function thispdf = eg_f(tau,mu,sigma,X)
    % PDF of RT(s) for a standard ex-Gaussian distribution with indicated parameters.
    rate = 1/tau;
    t1 = -X*rate + mu*rate + 0.5*(sigma*rate)^2;
    t2 = (X - mu - sigma^2*rate) / sigma;
    thispdf = rate*exp( t1 + log(normcdf(t2)) );
end
