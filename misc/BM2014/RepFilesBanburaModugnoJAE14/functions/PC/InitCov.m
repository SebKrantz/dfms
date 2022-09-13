function C = InitCov(T,P);
%___________________________________________________________________
%  PROCEDURE InitCov                                                
%  PURPOSE: Calculates the unconditional VCV of a stationary VAR(1) 
%  y = T y{1} + e,   e ~ (0,P)                                      
%                           P must be p.s.d.                                      
%                           T must have all evals < 1 in absval     
%  This is used for initializing the stationary part of the Kalman Filter                  
%  See Luetkepohl, p.22
%  OUTPUT   C with y ~ N(0,C) 
%___________________________________________________________________
% Check stationarity
  if max(abs(eig(T))) >= 1
     warning('InitCov: Matrix T has eigenvalues >= 1');
  end
 
% Check whether P is psd
  [V,D] = eig(P);
  if (min(diag(D))/max(diag(D)) < -1e-8) | (max(diag(D)) < 0)
     error('InitCov: Matrix P has eigenvalues < 0')
  end
 
% Calculate unconditional VCV
  m    = size(T,1);
  vecC = inv(eye(m^2) - kron(T,T)) * reshape(P,m^2,1);
  C    = reshape(vecC,m,m);
  C2   = C;
  C    = real(0.5*(C+C'));
    
   if abs(max(max(C2-C))) > 1e-8
      warning('InitCov: Large matrix adjustment ')
      disp('Max value of asymmetries')
      disp(num2str(abs(max(max(C2-C)))))
   end