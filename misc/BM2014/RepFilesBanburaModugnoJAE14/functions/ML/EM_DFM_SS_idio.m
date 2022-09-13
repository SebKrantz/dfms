function Res = EM_DFM_SS_idio(X,P)

% X: T x N panel of data
% P: structure containing settings of the model and estimation 

thresh = 1e-4;
r = P.r;%number of factors
p = P.p; %number of lags in the factor VAR
max_iter = P.max_iter; %maximal number of EM iterations
i_idio = P.i_idio; % logical vector Nx1 with 1's for variables with serially correlated idio

%--------------------------------------------------------------------------
% Preparation of the data
%--------------------------------------------------------------------------
[T,N] = size(X);

% Standardise x
Mx = nanmean(X);
Wx = (nanstd(X));
xNaN = (X-repmat(Mx,T,1))./repmat(Wx,T,1);
% xNaN = X;

%--------------------------------------------------------------------------
% Initial Conditions
%--------------------------------------------------------------------------

%Removing missing values (for initial estimators)

optNaN.method = 2; % Remove leading and closing zeros
optNaN.k = 3;

[A, C, Q, R, Z_0, V_0] = InitCond(xNaN,r,p,optNaN,i_idio);

% some auxiliary variables for the iterations
previous_loglik = -inf;
num_iter = 0;
LL = -inf;
converged = 0;

% y for the estimation is WITH missing data
y = xNaN';


%--------------------------------------------------------------------------
%THE EM LOOP
%--------------------------------------------------------------------------

%The model can be written as
%y = C*Z + e;
%Z = A*Z(-1) + v
%where y is NxT, Z is (pr)xT, etc

%remove the leading and ending nans for the estimation
optNaN.method = 3;
y_est = remNaNs_spline(xNaN,optNaN)';

while (num_iter < max_iter) & ~converged
    [C_new, R_new, A_new, Q_new, loglik, Z_0, V_0] = EMstep(y_est, A, C, Q, R, Z_0, V_0, r, p, i_idio);
    
    C = C_new;
    R(~i_idio,~i_idio) = R_new(~i_idio,~i_idio);
    A = A_new;
    Q = Q_new;

    % Checking convergence
    if num_iter>2
    [converged,decrease(num_iter+1)] = em_converged(loglik, previous_loglik, thresh,1);
    end
    
    LL = [LL loglik];
    previous_loglik = loglik;
    num_iter =  num_iter + 1;
end

%final run of the Kalman filter
Zsmooth = runKF(y, A, C, Q, R, Z_0, V_0)';
Res.x_sm = Zsmooth(2:end,:)*C';
Res.X_sm = repmat(Wx,T,1).*Res.x_sm+repmat(Mx,T,1);
Res.F = Zsmooth(2:end,:);

%--------------------------------------------------------------------------
%   Loading the structure with the results
%--------------------------------------------------------------------------
Res.C = C;
Res.R = R;
Res.A = A;
Res.Q = Q;
Res.Z_0 = Z_0;
Res.V_0 = V_0;
Res.r = r;
Res.p = p;
Res.Mx = Mx;
Res.Wx = Wx;

% Res.loglik   = LL;
% Res.num_iter = num_iter;
% Res.X_sm = X_sm;
% Res.converge = converged;
% Res.decrease = any(decrease);
% Res.decr_iter = min(find(decrease));
% decrease = any(decrease);

%--------------------------------------------------------------------------
%PROCEDURES
%--------------------------------------------------------------------------

function  [C_new, R_new, A_new, Q_new, loglik, Z_0, V_0] = EMstep(y, A, C, Q, R, Z_0, V_0, r, p, i_idio)

[n,T] = size(y);

% Compute the (expected) sufficient statistics for a single Kalman filter sequence.

%Running the Kalman filter with the current estimates of the parameters
[Zsmooth, Vsmooth, VVsmooth, loglik] = runKF(y, A, C, Q, R, Z_0, V_0);


EZZ = Zsmooth(1:r*p,2:end)*Zsmooth(1:r*p,2:end)'+sum(Vsmooth(1:r*p,1:r*p,2:end),3);                        %E(Z'Z)
EZZ_BB = Zsmooth(1:r*p,1:end-1)*Zsmooth(1:r*p,1:end-1)'+sum(Vsmooth(1:r*p,1:r*p,1:end-1),3); %E(Z(-1)'Z_(-1))
EZZ_FB = Zsmooth(1:r*p,2:end)*Zsmooth(1:r*p,1:end-1)'+sum(VVsmooth(1:r*p,1:r*p,:),3);%E(Z'Z_(-1)) 

EZZ2 = diag(diag(Zsmooth(r*p+1:end,2:end)*Zsmooth(r*p+1:end,2:end)'))+diag(diag(sum(Vsmooth(r*p+1:end,r*p+1:end,2:end),3)));                        %E(Z'Z)
EZZ_BB2 = diag(diag(Zsmooth(r*p+1:end,1:end-1)*Zsmooth(r*p+1:end,1:end-1)'))+diag(diag(sum(Vsmooth(r*p+1:end,r*p+1:end,1:end-1),3))); %E(Z(-1)'Z_(-1))
EZZ_FB2 = diag(diag(Zsmooth(r*p+1:end,2:end)*Zsmooth(r*p+1:end,1:end-1)'))+diag(diag(sum(VVsmooth(r*p+1:end,r*p+1:end,:),3)));%E(Z'Z_(-1)) 


A_new = A;
Q_new = Q;

A_new(1:r,1:r*p) = EZZ_FB(1:r,:) * inv(EZZ_BB);
Q_new(1:r,1:r) = (EZZ(1:r,1:r) - A_new(1:r,1:r*p)*EZZ_FB(1:r,:)') / T;

A_new(r*p+1:end,r*p+1:end) = EZZ_FB2 * diag(1./diag((EZZ_BB2)));
Q_new(r*p+1:end,r*p+1:end) = (EZZ2 - A_new(r*p+1:end,r*p+1:end)*EZZ_FB2') / T;

%E(Y'Y) & E(Y'Z) 
nanY = isnan(y);
y(nanY) = 0;

denom = zeros(n*r,n*r);
nom = zeros(n,r);
for t=1:T
    nanYt = diag(~nanY(:,t));
    denom = denom + kron(Zsmooth(1:r,t+1)*Zsmooth(1:r,t+1)'+Vsmooth(1:r,1:r,t+1),nanYt);
    nom = nom + y(:,t)*Zsmooth(1:r,t+1)'...
        -nanYt(:,i_idio)*(Zsmooth(r*p+1:end,t+1)*Zsmooth(1:r,t+1)'+Vsmooth(r*p+1:end,1:r,t+1));
    
end

vec_C = inv(denom)*nom(:);
C_new = C;
C_new(1:n,1:r) = reshape(vec_C,n,r);


R_new = zeros(n,n);
for t=1:T
    nanYt = diag(~nanY(:,t));
    R_new = R_new + (y(:,t)-nanYt*C_new*Zsmooth(:,t+1))*(y(:,t)-nanYt*C_new*Zsmooth(:,t+1))'...
        +nanYt*C_new*Vsmooth(:,:,t+1)*C_new'*nanYt...
        +(eye(n)-nanYt)*R*(eye(n)-nanYt);
end
R_new = R_new/T;
R_new = diag(diag(R_new));

% Initial conditions
V_0 = Vsmooth(1:r*p,1:r*p,1);
V2_0 = diag(diag(Vsmooth(r*p+1:end,r*p+1:end,1)));
V_0 = blkdiag(V_0,V2_0);

Z_0 = Zsmooth(:,1); 

%--------------------------------------------------------------------------

function [converged, decrease] = em_converged(loglik, previous_loglik, threshold, check_increased)
% EM_CONVERGED Has EM converged?
% [converged, decrease] = em_converged(loglik, previous_loglik, threshold)
%
% We have converged if the slope of the log-likelihood function falls below 'threshold', 
% i.e., |f(t) - f(t-1)| / avg < threshold,
% where avg = (|f(t)| + |f(t-1)|)/2 and f(t) is log lik at iteration t.
% 'threshold' defaults to 1e-4.
%
% This stopping criterion is from Numerical Recipes in C p423
%
% If we are doing MAP estimation (using priors), the likelihood can decrase,
% even though the mode of the posterior is increasing.

if nargin < 3, threshold = 1e-4; end
if nargin < 4, check_increased = 1; end

converged = 0;
decrease = 0;

if check_increased
    if loglik - previous_loglik < -1e-3 % allow for a little imprecision
        fprintf(1, '******likelihood decreased from %6.4f to %6.4f!\n', previous_loglik, loglik);
        decrease = 1;
    end
end

delta_loglik = abs(loglik - previous_loglik);
avg_loglik = (abs(loglik) + abs(previous_loglik) + eps)/2;
if (delta_loglik / avg_loglik) < threshold, converged = 1; end
 

%--------------------------------------------------------------------------

function [ A, C, Q, R, initZ, initV] = InitCond(xNaN,r,p,optNaN,i_idio)

[x, indNaN] = remNaNs_spline(xNaN,optNaN); 


OPTS.disp=0;

[T,N] = size(x);
% Eigenval decomp of cov(x) = VDV', only r largest evals
[ v, d ] = eigs(cov(x),r,'lm',OPTS);

% Observation equation   
C = [v zeros(N,r*(p-1))];

eyeN = eye(N);
eyeN(:,~i_idio) = [];
C=[C eyeN];

% Static predictions
f = x*v;

chi = x*v*v';
res=x-chi;
resNaN = res;
resNaN(indNaN) = nan;

R = diag(nanvar(resNaN));

%--------------------------------------------------------------------------
% Transition equation
%--------------------------------------------------------------------------

ii_idio = find(i_idio);
n_idio = length(ii_idio);
B = zeros(n_idio);
S = zeros(n_idio);

for i = 1:n_idio;
    R(ii_idio(i),ii_idio(i)) = 1e-04;

    res_i = resNaN(:,ii_idio(i));
    % number of leading zeros
    leadZero = max( find( (1:T)' == cumsum(isnan(res_i)) ) );
    endZero = max( find( (1:T)' == cumsum(isnan(res_i(end:-1:1))) ) );
    
    res_i = res(:,ii_idio(i));
    res_i(end-endZero:endZero) = [];
    res_i(1:leadZero) = [];
    
    B(i,i) = inv(res_i(1:end-1)'*res_i(1:end-1))*res_i(1:end-1)'*res_i(2:end,:);
    S(i,i) = cov(res_i(2:end)-res_i(1:end-1)*B(i,i));
end
 

% Estimate A & Q from stacked F(t) = A*F(t-1) + e(t);
z = f;
Z = [];
for kk = 1:p
    Z = [Z z(p-kk+1:end-kk,:)]; % stacked regressors (lagged SPC)
end;
z = z(p+1:end,:);
%% run the var chi(t) = A*chi(t-1) + e(t);

A = zeros(r*p,r*p)';
A_temp = inv(Z'*Z)*Z'*z;
A(1:r,1:r*p) = A_temp';
A(r+1:end,1:r*(p-1)) = eye(r*(p-1));

Q = zeros(p*r,p*r);
e = z  - Z*A_temp; % VAR residuals
Q(1:r,1:r) = cov(e); % VAR covariance matrix


% Initial conditions
initZ = zeros(size(Z,2)+n_idio,1); % %[randn(1,r*(nlag+1))]';
rp2 = (r*p)^2;
initV1 = reshape(inv(eye(rp2)-kron(A,A))*Q(:),r*p,r*p);
initV2 = diag(1./diag(eye(size(B,1))-B.^2)).*S;
initV = blkdiag(initV1, initV2);
A = blkdiag(A, B);
Q = blkdiag(Q, S);