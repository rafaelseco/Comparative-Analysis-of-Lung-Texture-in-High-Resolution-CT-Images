function data = GreyLevelDifferenceMethod(input,L)

    H = zeros(1,L+1);
   [lin,col] = size(input);
    
        for i=1:lin    % Lines.
            for j=1:(col-1) % Collumns.
                dif = abs(input(i,j)-input(i,j+1));  % Difference between the pixel and the next one.
                H(1,dif+1) = H(1,dif+1) + 1;
            end
        end
    
        h = H./sum(H);
    
        media=0;
    for k=0:L-1
        media = media + k*h(k+1);
    end
    
    variancia = 0;
    for k=0:L-1
        variancia = variancia + ((k - media)^2 * h(k+1));
    end
    
    ASMd = 0;
    IDMd = 0;
    CORd = 0;
    CONd = 0;
    
    for k=0:L-1
        ASMd = ASMd + (h(k+1))^2;
        IDMd = IDMd + (h(k+1)/(1+(k+1))^2);
        CORd = CORd + h(k+1)*((k+1)-media)/sqrt(variancia);
        CONd = CONd + h(k+1)*((k+1)^2);
    end
    
    data = [variancia ASMd IDMd CORd CONd];
end