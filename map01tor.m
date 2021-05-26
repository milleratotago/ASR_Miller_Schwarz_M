function r = map01tor(x01)
    % Convert any vector x of real numbers between eps and (1-eps) into a vector of real numbers r
    x01(x01<eps) = eps;
    x01(x01>1-eps) = 1-eps;
    r = x01 ./ (1-x01);
end

