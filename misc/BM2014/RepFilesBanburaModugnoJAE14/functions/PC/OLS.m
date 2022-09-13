function [b Sig sig_b t_val] = OLS(Y,Z,cflag,sflag)
%______________________________________________________________________________
% [b Sig sig_b t_val] = OLS(Y,Z,cflag,sflag)
% Estimates the model Y = Z*beta + u from Ordinary Least Squares
% Y might be nobs x k matrix: in this case k regressions are run
% Any missing data are skipped.
%
% INPUTS
%  Y          Series to be forecast                        (nobs x k)
%  Z          Regressors (estimates of factors)            (nobs x r)
%  cflag      = 1 -> Add constant to Z 
%  sflag      = 1 -> Standardise data before regression
%
% OUTPUT
%  b          Coefficients                                 (r x k)
%  Sig        Covmat of residuals - Choleski decomp!       (k x k)
%  sig_b      Standard error of estimates                  (r x k)
%  t_val      t-values                                     (r x k)
%______________________________________________________________________________
% Add constant & standardise 
  if  cflag
      Z = [ones(size(Z(:,1))) Z];
  end
  if  sflag
      Y = standardise(Y); 
      Z = standardise(Z);    
  end    

% Clear all obs with missing data
  ix    = ~isnan(sum(Y,2)) & ~isnan(sum(Z,2));
  Yc    =  Y(ix,:);
  Zc    =  Z(ix,:);
 
% Estimate beta & resids
  iZc   = inv(Zc'*Zc);
  b     = iZc * Zc' * Yc;
  U     = Y - Z*b;
  Sig   = cov(U(ix,:));
  
% Std dev and t-stat of b
  sig_b = sqrt(diag(iZc)*diag(Sig)');
  t_val = b ./ sig_b;
  
% Choleski decomp of Sig  
  Sig   = chol(Sig)';              % ! Adjust this later on !
 