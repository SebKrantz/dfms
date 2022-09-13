function [F V D IC] = Estim_SPC(X,r_max,IC_flag);
%__________________________________________________________________________
% function [F V D IC] = Estim_SPC(X,r_max,IC_flag);
% Estimate static principal components
% Possibly use Bai & Ng (2002) information criterion 
%
% X        Data matrix  (standardised & outlier corrected!)        [T x N]
% r_max    Max nr of factors 
% IC_flag  = 0   Use r_max for nr of factors (no estimation)
%          = 1   Use information criterion by Bai & Ng (2002)
%
% OUTPUT
% F                Estimated factors                             [T x r]
% V                Eigenvectors                                  [N x r]
% D                Eigenvalues                                   [r x r]
% IC               Value of criterion
%__________________________________________________________________________
% Check for nan
  if sum(sum(isnan(X))) > 0
     error('X contains missing values')
  end   
  
%__________________________________________________________________________
% Eigenval decomp of cov(X) = VDV', only r_max largest eigenvalues
  d.disp = 0;
  [V D]  = eigs(cov(X),r_max,'lm',d);	
  F      = X * V;                                 
  
%__________________________________________________________________________
% Select the number of factors
  if IC_flag == 0
     r     = r_max; 
     IC    = [];

  else    
   % Info crit penalty term 
     [T N] = size(X);
     PT    = log(min([N;T])) * (N+T)/(N*T);     

     IC    = zeros(1,r_max);
     for i = r_max:-1:1
         Xi       = X - F(:,1:i)*(V(:,1:i))';
         Sigma    = mean(sum(Xi.*Xi/T));         
         IC(i)    = log(Sigma) + i*PT;
     end
     
     [mv r]  = min(IC(1:r_max)');
     disp(['Nr of static factors = ' num2str(r)]);
   end
  
%__________________________________________________________________________
% Trim output 
  F  =  F(:,1:r);
  V  =  V(:,1:r);
  D  =  D(1:r,1:r);
  