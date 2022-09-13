function R = SEC(Ym,Model,Par,R,flag);
%______________________________________________________________________
% USAGE function R = SEC(Ym,Model,Par,R,flag);
% Signal extraction, i.e produces the signal for the obs vector based on
% either past or full sample information
% Always run after SKF and FIS
% The model is        y_t   = Z * a_t  + c1 + G * eps_t       
%                     a_t+1 = T * a_t  + c2 + H * eps_t       
% with HG'=0                                                         
%______________________________________________________________________
% INPUT  
%        Ym         Data                                 (nobs x n)  
%        Model      Name of function to build SSF        (string)
%        Par        Parameter structure to build SSF
%        R Estimates from Kalman filter FIS                                                          
%          R.AmU  : Estimates     a_t|t                    (nobs x m) 
%          R.PmU  : P_t|t   = Cov(a_t|t)               (nobs x m x m)       
%          R.AmT  : Estimates     a_t|T                    (nobs x m) 
%          R.PmT  : P_t|T   = Cov(a_t|T)               (nobs x m x m)
%        where m is the dim of state vector and t = 1 ...T is time
%        flag       signal based on estimates a_t|t (0) or a_t|T (1)
% OUTPUT 
%       R added
%         signal    signal for data based on (nobs x n)
%         sign_P    Standard deviation of signal
%______________________________________________________________________
% signal uses a function S = Model(Par,S,t), which builds the SSF
% (see SKF for more details)
%______________________________________________________________________
% Initalise
  eval(['S = ' Model '(Par,[],0);']);
  
  if size(S.Z,1) ~= size(Ym,2)
     error('Inappropriate length of data vector'); 
  end
  
% Dimensions  
  [nobs n]  = size(Ym);
  R.signal  = zeros(nobs,n);
  R.sign_P  = zeros(nobs,n);
  
% Go  
  for t = 1:nobs
      eval(['S = ' Model '(Par,S,' num2str(t) ');']);
            
      if flag
         R.signal(t,:)   = (S.Z * R.AmT(t,:)' + S.c1)';
         R.sign_P(t,:)   = sqrt(diag(S.Z * squeeze(R.PmT(t,:,:)) * S.Z'));
      else
         R.signal(t,:)   = (S.Z * R.AmU(t,:)' + S.c1)';  
         R.sign_P(t,:)   = sqrt(diag(S.Z * squeeze(R.PmU(t,:,:)) * S.Z'));
     end
  end
 