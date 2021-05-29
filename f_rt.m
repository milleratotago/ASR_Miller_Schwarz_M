function thispdf = f_rt(RT,Alpha,Beta,MuC,SigmaCab,Lambda,SigmaCba,SOA)
    % PDF(s) of overall RTs (con or inc) for a given set of parameter values.
    % Based on Formula A.6 in the manuscript
    % SigmaCab is sigma of stage C when a finishes before b, and SigmaCba with ba finishing order.
    % Lambda should be positive for inc and negative for con

    expAlphaSOA = exp(-Alpha*SOA);

    % q_s = Pr(A > B + SOA) = Pr(B finishes before A)
    q_s = Beta / (Alpha + Beta) * exp(-Alpha*SOA);

    % density when B finishes before A
    f_plus = eg_f(RT,1/(Alpha+Beta),MuC+Lambda,SigmaCba);

    % density when A finishes before B
    f_minus = (Alpha + Beta)   / (Alpha+Beta*(1-expAlphaSOA)) * eg_f(RT,1/Beta,MuC,SigmaCab) ...
            - Beta*expAlphaSOA / (Alpha+Beta*(1-expAlphaSOA)) * eg_f(RT,1/(Alpha+Beta),MuC,SigmaCab);

    thispdf = q_s * f_plus +  (1-q_s) * f_minus;

    % Due to numerical problems with some extreme parameter combinations, thispdf values
    % can sometimes be less than zero, which is of course impossible.
    % To circumvent that problem, we set negative pdf values to very small positive values.
    thispdf(thispdf<=0) = realmin;  % realmin is a MATLAB constant, with log(realmin) = -708.4
    
end
