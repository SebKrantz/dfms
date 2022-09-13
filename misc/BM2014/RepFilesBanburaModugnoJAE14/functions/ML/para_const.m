function Res = para_const(X,P, lag)


Z_0 = P.Z_0;
V_0 = P.V_0;
A = P.A;
C = P.C;
Q = P.Q;
R = P.R;
Mx = P.Mx;
Wx = P.Wx;

%--------------------------------------------------------------------------
% Preparation of the data
%--------------------------------------------------------------------------
[T,N] = size(X);
% Standardise x
xNaN = (X-repmat(Mx,T,1))./repmat(Wx,T,1);

y = xNaN';

%final run of the Kalman filter
[Zsmooth,P] = runKF_lag(y, A, C, Q, R, Z_0, V_0, lag);
Zsmooth=Zsmooth';
x_sm = Zsmooth(2:end,:)*C';
X_sm = repmat(Wx,T,1).*x_sm+repmat(Mx,T,1);

%--------------------------------------------------------------------------
%   Loading the structure with the results
%--------------------------------------------------------------------------
Res.P = P;
Res.X_sm = X_sm;

