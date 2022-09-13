function R = FIS(Ym,Model,Par,R);
%______________________________________________________________________
% USAGE function R = FIS(Ym,Model,Par,R)
% Fixed intervall smoother (see Harvey, 1989, p. 154)
% FIS returns the smoothed state vector AmT and its covar matrix PmT             
% The model is  
%                  y_t   = Z_t * a_t    + G_t * eps_t                   
%                  a_t+1 = T_t * a_t    + H_t * eps_t                   
% with HG'=0                                                         
% Use this in conjnuction with function SKF
%______________________________________________________________________
% INPUT  
%        Ym         Data                                 (nobs x n)  
%        Model      Name of function to build SSF        (string)
%        Par        Parameter structure to build SSF
%        R Estimates from Kalman filter SKF                                                          
%          R.Am   : Estimates     a_t|t-1                  (nobs x m) 
%          R.Pm   : P_t|t-1 = Cov(a_t|t-1)             (nobs x m x m)
%          R.AmU  : Estimates     a_t|t                    (nobs x m) 
%          R.PmU  : P_t|t   = Cov(a_t|t)               (nobs x m x m)       
% OUTPUT 
%        R Smoothed estimates added to above
%          R.AmT  : Estimates     a_t|T                    (nobs x m) 
%          R.PmT :  P_t|T   = Cov(a_t|T)               (nobs x m x m)
%        where m is the dim of state vector and t = 1 ...T is time
%______________________________________________________________________
% a) SKF uses a function S = Model(Par,S), which builds the state space
%    form. Par is a structure containing Parameter values and other
%    arbitrary information. Structure S contains matrices {Z, T, G, H}  
%    plus initial conditions A1 and P1 = cov(A1) for the state vector.
%    Use S=[] for the 1st call of Model and update S thereafter.
%    Model is passed to SKF as a text string
%______________________________________________________________________
  [nobs m]        = size(R.Am);
  R.AmT           = zeros(nobs,m);
  R.PmT           = zeros(nobs,m,m);
  R.AmT(nobs,:)   = squeeze(R.AmU(nobs,:))  ;
  R.PmT(nobs,:,:) = squeeze(R.PmU(nobs,:,:));
  
% Initalise & update to get S(nobs)
  eval(['S = ' Model '(Par,[],0);']);
  eval(['S = ' Model '(Par,S,' num2str(nobs) ');']);
  
  for t = nobs-1:-1:1 
      PmU = squeeze(R.PmU(t,:,:));
      Pm1 = squeeze(R.Pm(t+1,:,:));
      P_T = squeeze(R.PmT(t+1,:,:));
      
      eval(['S = ' Model '(Par,S,' num2str(t) ');']);
      C  = PmU * S.T' * pinv(Pm1);
   
      R.AmT(t,:)   = R.AmU(t,:) + (C * (R.AmT(t+1,:) - R.Am(t+1,:))')' ; 
      R.PmT(t,:,:) = PmU        +  C * (P_T - Pm1) * C'; 
      R.Pstar(t,:,:)=C;
  end
 