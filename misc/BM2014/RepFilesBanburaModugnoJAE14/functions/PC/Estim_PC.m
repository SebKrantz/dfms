function P = Estim_PC(F,P,n)
%__________________________________________________________________________
% function P = Estim_PC(F,P,n)
% Estimates the dynamic (2nd) equation of the factor model
%
%         X_t = C F_t            +   xi_t,       xi_t ~ WN(0,R)
%         F_t = Sum_p(A_i F_{t-i}) +  u_t,        u_t ~ WN(0,Q)
%
% which is a VAR of order P.p and with rank(Q) = P.q
%
% INPUT
%  F             Estimates of factors                          [nobs x r]
%  P.crit           = 1 Use information criteria for p & q
%                   = 0 Take values P.p & P.q from P
%  P.q           Nr of dyn factors                          
%  P.p           Nr of lags in factor VAR   
%  n             Nr of series in original data (used for static PC)
%
% OUTPUT
%   P            Extended with parameter matrices A Q H &
%                  possibly adjusted values of q & p (qflag > 0)                  
%__________________________________________________________________________
  [nobs r] = size(F);

%__________________________________________________________________________
% Find optimal lag using function V_AR
     p         = P.p; 
%__________________________________________________________________________
% Stacked AR regression F(t) = A*F(t-1) + e(t)
  I  = eye(r*p);    
  A  = [   zeros(r,r*p)  ; ... 
        I(1:end-r,1:end) ];    
      
% Form expl vars Z for AR regression
% Estimate stacked F(t) = A*F(t-1) + e(t)  
  Z = [];
  for i = 1:p
      Z = [Z F(p-i+1:end-i,:)];
  end
  z            = F(p+1:end,:);             
  A1           = inv(Z'*Z)*Z'*z;          
  E            = z  - Z*A1; 
  A(1:r,1:r*p) = A1';
  
%__________________________________________________________________________
% Find optimal q as from Bai & Ng (2005)
      q = min([P.q r]);
  
%__________________________________________________________________________
% Form reduced rank covariances Q = H*H' from resids E
  if r > 1 
     d.disp  = 0;
     [D M]   = eigs(cov(E),q,'lm',d);
     D       = D*diag(sign(D(1,:)));
  else
     D = 1;
     M = cov(E);
  end
  H          = zeros(p*r,p*r);
  Q          = zeros(p*r,p*r);
  H(1:r,1:q) = D*sqrt(M);
  Q(1:r,1:r) = D*M*D';

%__________________________________________________________________________
  P.A  = A; 
  P.H  = H;
  P.Q  = Q;
  P.p  = p;
  P.q  = q;
  