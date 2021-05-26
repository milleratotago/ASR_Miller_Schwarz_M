function x01 = maprto01(realx)
    % Convert any vector of real numbers x into a vector of numbers between eps and (1-eps)
    absx = abs(realx);
    x01 = absx ./ (1+absx);
    x01(x01<eps) = eps;
    x01(x01>1-eps) = 1-eps;
end

