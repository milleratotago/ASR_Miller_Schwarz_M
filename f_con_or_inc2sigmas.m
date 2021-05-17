function thispdf = f_con_or_inc2sigmas(TauA,TauB,MuC,SigmaCab,Lambda,SigmaCba,RT,varargin)
    % PDF(s) of RT_con(s) or RT_inc(s) for a given set of parameter values.
    % SigmaCab is sigma of stage C when a finishes before b, and SigmaCba with ba finishing order.
    % Lambda should be positive for inc and negative for con
    % Optional parameter in varargin is SOA; assume 0 if not used.

    if numel(varargin) == 0
        SOA = 0;
    else
        SOA = varargin{1};
    end

    Alpha = 1/TauA;
    Beta = 1/TauB;
    expAlphaSOA = exp(-Alpha*SOA);

    % q_s = Pr(A > B + SOA) = Pr(B finishes before A)
    q_s = Beta / (Alpha + Beta) * exp(-Alpha*SOA);

    % density when B finishes before A
    f_plus = eg_f(1/(Alpha+Beta),MuC+Lambda,SigmaCba,RT);

    % density when B finishes after A
    f_minus = (Alpha + Beta)   / (Alpha+Beta*(1-expAlphaSOA)) * eg_f(1/Beta,MuC,SigmaCab,RT) ...
            - Beta*expAlphaSOA / (Alpha+Beta*(1-expAlphaSOA)) * eg_f(1/(Alpha+Beta),MuC,SigmaCab,RT);

    thispdf = q_s * f_plus +  (1-q_s) * f_minus;

    % Due to numerical problems with some extreme parameter combinations, thispdf values
    % can sometimes be less than zero, which is of course impossible.
    % To circumvent that problem, we set negative pdf values to very small positive values.
    thispdf(thispdf<=0) = realmin;  % realmin is a MATLAB constant, with log(realmin) = -708.4
    
end
