function [quanCon, quanInc, quanMn, quanDelta, quanConPDF, quanIncPDF] = ASRplots(RTcon,RTinc,cumProbs)
    % Plot PDFs, CDFs, and deltas for the indicated congruent and incongruent RT distributions.
    % Evaluate the PDFs, CDFs, and deltas at quantile points indicated by the vector in cumProbs.
    
    %% Preliminary computations:
    
    [quanCon, quanInc, quanMn, quanDelta] = ASRquans(RTcon,RTinc,cumProbs);
    
    quanConPDF = zeros(size(cumProbs));
    quanIncPDF = quanConPDF;
    
    for i = 1:numel(cumProbs)
        quanConPDF(i) = RTcon.PDF(quanCon(i));
        quanIncPDF(i) = RTinc.PDF(quanInc(i));
    end
    
    %% Figure
    
    subplot(3,1,1);
    plot(quanCon,quanConPDF);
    hold on
    plot(quanInc,quanIncPDF);
    legend('Congruent','Incongruent');
    xlabel('RT')
    ylabel('PDF');
    
    subplot(3,1,2);
    plot(quanCon,cumProbs);
    hold on
    plot(quanInc,cumProbs);
    legend('Congruent','Incongruent');
    xlabel('RT')
    ylabel('CDF');
    
    subplot(3,1,3);
    plot(quanMn,quanDelta);
    xlabel('RT')
    ylabel('\Delta');
    
end

