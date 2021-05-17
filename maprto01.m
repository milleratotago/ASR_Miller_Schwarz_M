function x01 = maprto01(realx)
    % Convert any real number x into a number between 0-1
    absx = abs(realx);
    x01 = absx / (1+absx);
end

