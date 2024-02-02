% @veronica
% Performs the quantization of the matrix I, with NL  gray levels. The
% maximum value of I is G(2) and the minimum is G(1). 

function SI=scaleNL (I, NL, GL)
I=double(I);
if GL(2) == GL(1)
  SI = ones(size(I))-1;
else
  slope = (NL-1) / double(GL(2) - GL(1));
  intercept = 1 - double((slope*(GL(1))));
  SI = round(imlincomb(slope,I,intercept,'double'))-1;
end
end