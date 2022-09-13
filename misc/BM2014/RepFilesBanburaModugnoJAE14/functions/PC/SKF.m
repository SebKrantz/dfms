function R = SKF(Ym,Model,Par)
%______________________________________________________________________
% USAGE R = SKF(Ym,Model,Par)
% Kalman filter for stationary systems with time-varying system matrices
% and missing data.
% Version for block correlations in matrix G*G'
%
% The model is        y_t   = Z(t) * a_t  + c1(t) + G(t) * eps_t       
%                     a_t+1 = T(t) * a_t  + c2(t) + H(t) * eps_t       
%
% with H'G = 0  & diagonal matrix G G'                                                         
%
%                            I M P O R T A N T
%  A special expression is used to calculate the inverse of F_t, as 
%  the dimension of the obs eq is high. This requires S.G*S.G' to be  
%  of a certain structure. The first S.nd rows/cols are assumed to
%  be diagonal, whereas now assumptions are made on the second part.

%  If this is not the case, change the code below & compute 
%  inv(F_t) as iF = inv(Z*P*Z' + G*G')
%______________________________________________________________________
% INPUT  
%        Ym         Data                                 (nobs x n)  
%        Model      Name of function to build SSF        (string)
%        Par        Parameter structure to build SSF
% OUTPUT 
%        R.Am       Predicted state vector  A_t|t-1      (nobs x m)  
%        R.AmU      Filtered  state vector  A_t|t        (nobs x m)  
%        R.Pm       Predicted covariance of A_t|t-1      (nobs x m x m)  
%        R.PmU      Filtered  covariance of A_t|t        (nobs x m x m)  
%        R.Q        Sum of v'F^{-1}v                                   
%        R.SF       v' inv(F) v
%        R.loglik   Value of likelihood function
%______________________________________________________________________
% a) SKF uses a function S = Model(Par,S), which builds the state space
%    form. Par is a structure containing Parameter values and other
%    arbitrary information. Structure S contains matrices {Z, T, G, H}  
%    plus initial conditions A1 and P1 = cov(A1) for the state vector.
%    Use S=[] for the 1st call of Model and update S thereafter.
%    Model is passed to SKF as a text string
% b) Function MissData checks for missing series in obs i (Ym(i,:))
%    it reduces the observation equation accordingly. If observation i
%    is empty, a pure forecasting step is done. Otherwise the KF  
%    recursion is applied to the reduced obs eq.
%______________________________________________________________________
% Initalise
  eval(['S = ' Model '(Par,[],0);'])
  if size(S.Z,1) ~= size(Ym,2)
      disp(['Data vector : ' num2str(size(Ym,2))]);
      disp(['Z           : ' num2str(size(S.Z,1))]);
      error('Data vector does not fit Z'); 
  end
  
% Output structure & dimensions
  [n m] = size(S.Z);
  nobs  = size(Ym,1);
  
  R.sF  = 0;
  R.Q   = 0;
  R.Am  = nan*zeros(nobs,m);   R.Pm  = nan*zeros(nobs,m,m);
  R.AmU = nan*zeros(nobs,m);   R.PmU = nan*zeros(nobs,m,m);
   
%______________________________________________________________________   
  A  = S.A1;
  P  = S.P1;
  Au = zeros(size(A));
  
  for t = 1:nobs
      R.Am(t,:)   = A;
      R.Pm(t,:,:) = P;

    % Obtain SSF and D = inv(GG') 
      eval(['S = ' Model '(Par,S,' num2str(t) ');']);
      [y,Z,G,c1,L] = MissData(Ym(t,:)',S);

      if isempty(y)
         Au = A;
         Pu = P;
         A  = S.T*A      + S.c2;
         P  = S.T*P*S.T' + S.H*S.H';  
         iF = zeros(n,n);
         K  = zeros(m,n);         
      else
          
       % Compute inv(F_t)
         PZ = P*Z';
         GG = G*G';

         if sum(sum(GG-diag(diag(GG)))) > 0 | sum(find(diag(GG)==0)) > 0
            iF  = inv(Z*PZ + GG);
         else
            D   = diag(1./diag(GG));
            iF  = D - D*Z*inv(eye(m)+PZ*D*Z)*PZ*D;
         end
         if sum(sum(abs(iF*(Z*P*Z'+G*G') - eye(size(iF))))) > n^2*1e-6
            error('Inversion of matrix F in did not work')
            return
         end
       
       % Kalman gain K_t  
         PZF = PZ*iF;
         K   = S.T*PZF;
         
       % Au = A_t|t   & Pu = P_t|t
         V   = y - Z*A  -  c1;
         Au  = A  + (PZF * V);           
         Pu  = P  - (PZF * PZ');
             
       % A  = A_t+1|t & P  = P_t+1|t  
         A   =  S.T*A + K*V + S.c2;
         P   = (S.T-K*Z)*P*S.T' + S.H*S.H';
         P   =  0.5 * (P+P');           
       
       % Likelihood  
         R.Q  = R.Q  + V'*iF*V;
         R.sF = R.sF - log(det(iF)); 
         
       % Restore structure of iF and K  
         iF   = L*iF*L';
          K   = K*L';
     end
     
     R.AmU(t,:)    = Au;  
     R.PmU(t,:,:)  = Pu;
      R.iF(t,:,:)  = iF;
       R.K(t,:,:)  = K;
  end % t 
  
% Likelihood  
  R.loglik = 0.5 * (R.sF + R.Q);
