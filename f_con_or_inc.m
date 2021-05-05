function thispdf = f_con_or_inc(TauA,TauB,MuC,SigmaC,Lambda,RT,varargin)
    % PDF(s) of RT_con(s) or RT_inc(s) for a given set of parameter values.
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
    f_plus = eg_f(1/(Alpha+Beta),MuC+Lambda,SigmaC,RT);

    % density when B finishes after A
    f_minus = (Alpha + Beta)   / (Alpha+Beta*(1-expAlphaSOA)) * eg_f(1/Beta,MuC,SigmaC,RT) ...
            - Beta*expAlphaSOA / (Alpha+Beta*(1-expAlphaSOA)) * eg_f(1/(Alpha+Beta),MuC,SigmaC,RT);

    thispdf = q_s * f_plus +  (1-q_s) * f_minus;
    
end

