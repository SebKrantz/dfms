function [R_KF] = RUN_DFM_ECB(X,y,Date,Q)

n = size(X,2);

%__________________________________________________________________________
% Static principal components
if Q.outl
    Xe        = Outliers(z01(X(1:end-Q.cutE,:)),2);
else
    Xe        = Outliers(z01(X(1:end-Q.cutE,:)),5);
end

[F V]      = Estim_SPC(Xe,max(Q.r),0);
FQ        = filter([1/3 1/3 1/3],1,F);
FQ(1:2,:) = nan;


Xf = z01(X);
if Q.outl
    Xf          = Outliers(Xf,0);
end
yf  = y;

Q.Date      = Date;
Tf          = size(Xf,1);


% Estimate dynamic model
Fj    = F(:,1:Q.r);
Vj    = V(:,1:Q.r);
Q     = Estim_PC(Fj,Q,n);
Q.C   = [Vj zeros(n,Q.r*(Q.p-1))];
Q.R   = diag(diag(cov(Xe-Fj*Vj')));
[Q.beta Q.s] = OLS(y(3:end-Q.cutE,:),FQ(3:end,1:Q.r),1,0);

% Correct for a too small variance of the idiosyncratic component
% Otherwise the KF would break down
if min(diag(Q.R)) < 1e-2;
    LL         = find(diag(Q.R)<1e-2);
    for i = 1:length(LL)
        Q.R(LL(i),LL(i)) = 1e-2;
    end
end

% Kalman filter & smoother
R_KF  = SKF([Xf yf],Q.Model,Q);
R_KF  = FIS([Xf yf],Q.Model,Q,R_KF);
R_KF  = SEC([Xf yf],Q.Model,Q,R_KF,1);


