function data = featuresFirstOrder(ROI,L)

H = zeros(1,L);
[alt,larg] = size(ROI);
for I = 1:alt
  for J = 1:larg
       index=ROI(I,J)+1; 
       H(index)=H(index)+1;
  end
end


h=H./sum(H);

media = 0;

for k=0:L-1
    media = media + k*h(k+1);
end

variancia = 0;
for k=0:L-1
    variancia = variancia + ((k - media)^2 * h(k+1));
end

simetria = 0;
for i=0:L-1
    simetria = simetria + ((i - media)^3 * h(i+1));
end
simetria = simetria/(sqrt(variancia)^3); 


curtose = 0;
for i=0:L-1
    curtose = curtose + (h(i+1)*(i-media)^4);
end
curtose = curtose/(sqrt(variancia)^4) - 3;

energia =0;
for k=0:L-1
    energia = energia + (h(k+1))^2;
end

entropia = 0;
for k=0:L-1
    if(h(k+1)~=0)
        entropia = entropia + h(k+1)*log2(h(k+1));
    end
end
entropia = - entropia;
  
data = [media,variancia,simetria,curtose,energia,entropia];
end