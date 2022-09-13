function [Y, A, const, Sig, IC] = V_AR(Y,P,h)
%______________________________________________________________________________
% [Y, A, const, Sig, AIC] = V_AR(Y,lags,h,cflag)
% Estimates an vector autoregressive model Y = c + A(L)*Y(-1) + u from OLS
%
% INPUTS
%  Y               Series                                        (nobs x n)
%  P.cflag      = 1 -> Add constant to VAR 
%
%  P.IC         Info crit for lag length selection
%                  'AIC'  Akaike
%                  'SIC'  Schwartz
%                  'FIX'  fixed lag length
%
%  P.lags          if IC = 'AIC', 'SIC' Maximum nr of lags in VAR                             (integer)
%                  if IC = 'FIX'        Nr of lags in VAR
%
%  h           Nr of steps in forecasting                        (integer)
%
% OUTPUT
%  Y          Series extended with dynamic forecast h steps     (nobs+h x n)
%  A          Coefficient matrices                              (n x n x lags)
%  const      Estimated constant                                (n x 1)
%  Sig        Covmat of residuals                               (n x n)
%  IC         Akaike / Schwartz information criterion
%______________________________________________________________________________

% Checks
  if sum(strcmp(P.IC,upper({'AIC','SIC','FIX'}))) == 0 
     error('Input P.IC is invalid')
  end
  if P.lags <= 0
     error('Lags must be positive integer')
  end    
      
% Estimate
  if strcmp(P.IC,upper('FIX'))
          lags  = P.lags;
         [b,IC] = Estimate_VAR(Y,P.lags+1,lags,P.cflag,'');
  else
     IC_min = inf;
     for i = 1:P.lags
         [b,IC] = Estimate_VAR(Y,P.lags+1,i,P.cflag,P.IC);
         if IC < IC_min
            lags     = i;
            IC_min  = IC;
         end
     end
     [b, IC, A, const, Sig] = Estimate_VAR(Y,P.lags+1,lags,P.cflag,P.IC);
  end
  
% Dynamic Forecasts
  [nobs n]  = size(Y);

  for s = 1:h
      Z_nxt = [];
      for jp = 0:lags-1
          Z_nxt = [Z_nxt Y(end-jp,:)];
      end
          
      if  P.cflag
          Z_nxt = [Z_nxt 1];
      end
      
      Y_nxt  = Z_nxt * b;
      Y = cat(1,Y,Y_nxt);
  end    
  
  
 function [b, IC, A, const, Sig] = Estimate_VAR(Y,cut,lags,cflag,Iflag)
%______________________________________________________________________________
% Estimates the above VAR for a fixed lag number
% INPUTS
%  Y               Series                                        (nobs x n)
%  cut             cut first cut obs in data
%  lags            Nr of lags in VAR
%  cflag           = 1 -> Add constant to VAR 
%  Iflag           Info crit for lag length selection
%______________________________________________________________________________

% Create lags & add constant
  Z = [];  
  for kk = 1:lags
      Z = [Z Y(lags-kk+1:end-kk,:)]; % stacked regressors (lagged SPC)
  end;
  if  cflag
      Z = [Z ones(size(Z,1),1)];
  end
  Z = [nan(lags,size(Z,2)); Z];
  Z = Z(cut+1:end,:);
  Y = Y(cut+1:end,:);
  
% Clear all obs with missing data
  ix    = ~isnan(sum([Y Z],2));
  Yc    =  Y(ix,:);
  Zc    =  Z(ix,:);

% Estimate beta & resids
  iZc   = inv(Zc'*Zc);
  b     = iZc * Zc' * Yc;
  U     = Y - Z*b;
  Sig   = cov(U(ix,:));
  
% Calculate AIC / SIC
  [nobs n] = size(Yc);
   if     strcmp(Iflag,upper('AIC'))
          IC  = log(det(Sig)) + (n^2*lags + n*cflag) * (2/nobs); 
   elseif strcmp(Iflag,upper('SIC'))    
          IC  = log(det(Sig)) + (n^2*lags + n*cflag) * (log(nobs)/nobs); 
   else
          IC  = nan; 
   end    
          
% Coefficent matrices  
  A  = nan(n,n,lags);
  for j = 1:lags
      A(:,:,j) = b((j-1)*n+1:j*n,:);
   end
 
  if cflag 
     const = b((lags+1)*n,:)';
  else
      const = [];
  end    