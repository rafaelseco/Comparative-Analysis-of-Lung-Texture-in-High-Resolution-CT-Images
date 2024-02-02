function result = entropyProject(Matrix,L)
%ENTROPY Summary of this function goes here
%   Detailed explanation goes here

matrix = Matrix./sum(Matrix(:));

Entropy=0;
for lin=1:L
    for col=1:L
        if(matrix(lin,col)~=0)
            Entropy = Entropy + matrix(lin,col)*log2(matrix(lin,col));
        end
    end
end

result= -Entropy;
end

