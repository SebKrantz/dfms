function [y0, meanv, stdev] = z01(y)
%______________________________________________________________________________
% function [y0, m, S] = z01(x)
% Produces standardised data set x0 with mean 0 and variance 1 per series
% Works also with missing data
%
% INPUT
%     Data               y   [nobs x n]
% OUTPUT
%     Standardised data  y0   [nobs x n]
%     Means                      [1 x n]
%     Standard deviations        [1 x n]
%______________________________________________________________________________
 [nobs,n]  = size(y);
  meanv    = zeros(1,n);
  stdev    = zeros(1,n);  
 
 for k = 1:n
      e        = ~isnan(y(:,k));
      if any(e)
          stdev(k) =    std(y(e,k));	 	           
          meanv(k) =  mean(y(e,k));
      else
          stdev(k) = nan; 	           
          meanv(k) = nan;
      end
 end
 
 y0 = (y - ones(nobs,1)*meanv) ./ (ones(nobs,1)* stdev);		