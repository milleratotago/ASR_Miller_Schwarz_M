function thispdf = eg_f(tau,mu,sigma,X)
    % PDF of RT(s) for a standard ex-Gaussian distribution with indicated parameters.
    rate = 1/tau;
    t1 = -X*rate + mu*rate + 0.5*(sigma*rate)^2;
    t2 = (X - mu - sigma^2*rate) / sigma;
    thispdf = rate*exp( t1 + log(normcdf(t2)) );
    % The above is the theoretical definition of the ex-Gaussian pdf,
    % but there are numerical problems if t1 > 708 or so, because then exp(t1) overflows.
    % That happens when tau is small relative to sigma (so (sigma*rate)^2 is large).
    % In that case, though, the exG density is very close to a normal with mean mu+tau & variance sigma^2+tau^2,
    % so we can just use that corresponding normal pdf in those cases to avoid numerical problems:
    t1bad = t1 > 708;  % vector indicating too-large t1's for which we should use the normal approximation.
    thispdf(t1bad) = normpdf(X(t1bad),mu+tau,sqrt(sigma^2+tau^2));
end
