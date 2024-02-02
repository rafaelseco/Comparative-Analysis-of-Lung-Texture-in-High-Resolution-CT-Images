function h = Histogram(roiGS,L,k,graphtitle)

%   Cálculo do Histograma e do Normalizado.

    [alt, larg] = size(roiGS);
    H= zeros(1,L);
    for I = 1:alt
        for J = 1:larg
            index=roiGS(I,J)+1; %Para estar contido em [1-256]
            H(index)=H(index)+1;
        end
    end

%      figure(k)
%         
%           h = H./sum(H);
%           h = h(1:1:L);
%           horz = 1:1:L;
%           bar(horz, h), 
%           axis([0 L 0 max(h)])
%           xlabel('Grey Levels');
%           ylabel('Number of occurrences');
%           title(graphtitle)
         


end