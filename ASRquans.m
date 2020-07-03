function [quanCon, quanInc, quanMn, quanDelta] = ASRquans(RTcon,RTinc,cumProbs)
    % Find the quantile points indicated by the vector in cumProbs.
    %   for the indicated congruent and incongruent RT distributions.
    % Evaluate the PDFs, CDFs, and deltas at these quantiles.
    
    quanCon = zeros(size(cumProbs));
    quanInc = quanCon;
    
    for i = 1:numel(cumProbs)
        quanCon(i) = RTcon.InverseCDF(cumProbs(i));
        quanInc(i) = RTinc.InverseCDF(cumProbs(i));
    end
    
    quanMn = (quanCon + quanInc) / 2;
    quanDelta = quanInc - quanCon;
    
end

