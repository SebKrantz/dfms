function [F_hat,F_pc,F_kal,num_iter] = dynFA(X,q,r,p,max_iter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "A Quasi?Maximum Likelihood Approach for Large, Approximate Dynamic Factor Models," 
% The Review of Economics and Statistics, MIT Press, vol. 94(4), pages 1014-1024, November 2012.
% Catherine Doz, Universite' Cergy-Pontoise
% Domenico Giannone, Universite' Libre de Bruxelles, ECARES and CEPR
% Lucrezia Reichlin, London Business School and CEPR 
%
%
% Programs are also available at: http://homepages.ulb.ac.be/~dgiannon/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% DynFA:             extracts the unobservable factors using three different methods 
%                         
%       - QML:      Max Likelihood estimates using the Expectation Maximization (EM) algorithm 
%                    (Doz, Giannone and Reichlin, 2012) 
%                         
%       - TWO STEP: Principal components + Kalman filtering 
%                   Doz, Catherine & Giannone, Domenico & Reichlin, Lucrezia, 2011.
%                   "A two-step estimator for large approximate dynamic factor models based on Kalman filtering," 
%                   Journal of Econometrics, Elsevier, vol. 164(1), pages 188-205, September.
% 
%       - PC:       principal components 
%                         
% INPUTS
% X - matrix of observable variables
% r - # of static factors
% q - # of dynamic factors
% p - # length of ar filter on common factors
% max_iter - max # of iterations in estimation
%
% OUTPUTS
% F_hat -   factors from QML
% F_pc  -   factors using principal components
% F_kal -   factors using from two steps


thresh = 1e-4;
OPTS.disp = 0;
[T,N] = size(X);

Mx = mean(X);
Wx = (std(X));
% x is the matrix of standardized observable variables
x = (X-kron(ones(T,1),Mx))*diag(1./Wx);

% the number of static factors cannot be great of the dynamic ones
if r < q;
    error('q has to be less or equal to r')
end

nlag = p-1;
                                                                                                    % nlag = p-1;

A_temp = zeros(r,r*(nlag + 1))';
I = eye(r*(nlag+1),r*(nlag+1));
A = [A_temp';I(1:end-r,1:end)];

Q = zeros((nlag+1)*r,(nlag+1)*r);
Q(1:r,1:r) = eye(r);

OPTS.disp=0;

%extract the first r eigenvectors and eigenvalues from cov(x)
[ v, d ] = eigs(cov(x),r,'lm',OPTS);

chi = x*v*v';                       % common component

d = eye(r);

F = x*v;

F_pc = F;                           % factors using principal components

if p > 0    
    z = F;
    Z = [];
    for kk = 1:p
        Z = [Z z(p-kk+1:end-kk,:)]; % stacked regressors (lagged SPC)
    end;
    z = z(p+1:end,:);
    % run the var chi(t) = A*chi(t-1) + e(t);
    A_temp = inv(Z'*Z)*Z'*z;        % OLS estimator of the VAR transition matrix
    A(1:r,1:r*p) = A_temp';
    e = z  - Z*A_temp;              % VAR residuals
    H = cov(e);                     % VAR covariance matrix
    if r > q    % if s is different from 0 we have q dynamic factors
        
        %extract the first q eigenvectors and eigenvalues from cov(e)
        [ P , M ] = eigs(H,q,'lm',OPTS);    
        P = P*diag(sign(P(1,:)));
        u_orth = e*P*(M^-.5);       % extract the common shocks
        e_pc = e*P*P';
        Q(1:r,1:r) = P*M*P';        % variance of the VAR shock when s>0
    else        
        Q(1:r,1:r) = H;             % variance of the VAR shock when s=0
    end;
end;

R = diag(diag(cov(x-chi)));         % R diagonal

z = F;
Z = [];
for kk = 0:nlag
    Z = [Z z(nlag-kk+1:end-kk,:)];  % stacked regressors (lagged SPC)
end;

initx = Z(1,:)';                    %initial state mean                                                            


% initial state covariance
initV = reshape(pinv(eye((r*(nlag+1))^2)-kron(A,A))*Q(:),r*(nlag+1),r*(nlag+1)); % eye(r*(nlag+1));


C = [v zeros(N,r*(nlag))];


% initialize the estimation and define ueful quantities

previous_loglik = -inf;

loglik = 0;
num_iter = 0;
LL = [];

os = size(C,1);     % number of cross sections ( N )
ss = size(A,1);     % number of factors ( r )
y = x';

LL = -inf;

converged = 0;

% estimation of the factors with the Kalman filter using as initial values
% for A, C, Q, R, initx, initV the ones computed with the principal
% components
[xitt,xittm,Ptt,Pttm,loglik_t]=K_filter(initx,initV,x,A,C,R,Q);
[xsmooth, Vsmooth, VVsmooth]=K_smoother(A,xitt,xittm,Ptt,Pttm,C,R);

F_kal =  xsmooth(1:r,:)';


% factors estimation with the EM algorithm

%repeat the algorithm until convergence
while (num_iter < max_iter) & ~converged
    
    %%% E step : compute the sufficient statistics 
    
    % In our model the sufficient statistics are
    % delta = sum_t=1^T (x_t * f'_t)
    % gamma = sum_t=1^T (f_t * f'_t)   
    % gamma1 = sum_t=2^T (f_t-1 * f'_t-1)
    % gamma2 = sum_t=2^T (f_t * f'_t)
    % beta = sum_t=1^T (f_t * f'_t-1)
    % P1sum  variance of the initial state
    % x1sum  expected value of the initial state 
    
    
    num_iter;
    
    % initialize to zero all the sufficient statistics
    delta = zeros(os, ss);      
    gamma = zeros(ss, ss);         
    gamma1 = zeros(ss, ss);     
    gamma2 = zeros(ss, ss);    
    beta = zeros(ss, ss);      
    P1sum = zeros(ss, ss);     
    x1sum = zeros(ss, 1);      
    loglik = 0;
    
    % use the function Estep to update the expected sufficient statistics
    % note that at the first iteration  we use as initial values for A, C, Q, 
    % R, initx, initV the ones computed with the principal components
    [beta_t, gamma_t, delta_t, gamma1_t, gamma2_t, x1, V1, loglik_t,xsmooth] = ...
        Estep(y, A, C, Q, R, initx, initV);
    
    % fix the expected sufficient statistics equal to the one computed with
    % the function Estep
    beta = beta + beta_t;
    gamma = gamma + gamma_t;
    delta = delta + delta_t;
    gamma1 = gamma1 + gamma1_t;
    gamma2 = gamma2 + gamma2_t;
    P1sum = P1sum + V1 + x1*x1';    
    x1sum = x1sum + x1;
    
                                                                            
    % update the loglikelihood                                                                           
    loglik = loglik + loglik_t;
    
    % update the counter for the iterations
    num_iter =  num_iter + 1;
    
    
    %%% M step 
    % compute the parameters of the model as a function of the sufficient
    % statistics (computed with the function Estep)
    
    % The formulas for the parameters derive from the maximum likelihood
    % method. In the EM algorithm we substitute in the ML estimator the 
    % sufficient statistics (computed in the E step) and then iterate the 
    % procedure till the maximization of the likelihood
    
    % C = (sum_t=1^T x_t*f'_t)* (sum_t=1^T f_t*f'_t)^-1 
    % substituting for the sufficient statistics
    
    C(:,1:r) = delta(:,1:r) * pinv(gamma(1:r,1:r));
    
    if p > 0
    
        % A = (sum_t=2^T f_t*f'_t-1)* (sum_2=1^T f_t-1*f'_t-1)^-1 
        Atemp = beta(1:r,1:r*p) * inv(gamma1(1:r*p,1:r*p));
        A(1:r,1:r*p) = Atemp;
        
        % Q = ( (sum_t=2^T f_t*f'_t) - A * (sum_2=1^T f_t-1*f'_t) )/(T-1) 
        H = (gamma2(1:r,1:r) - Atemp*beta(1:r,1:r*p)') / (T-1);
        if r > q
            [ P , M ] = eigs(H,q,'lm',OPTS);
            P = P*diag(sign(P(1,:)));
            u_orth = e*P*(M^-.5);       % extracting the common shocks
            e_pc = e*P*P';
            Q(1:r,1:r) = P*M*P';        % Q if s>0
        else 
            Q(1:r,1:r) = H;             % Q if s=0
        end;
    end;
                                                                                                        % eigs(Q)
    
    % R = ( sum_t=1^T (x_t*x'_t) - C * f_t*x'_t) )/T 
    R = (x'*x - C*delta')/T;
    
    RR = diag(R); RR(RR<1e-7) = 1e-7; R = diag(RR); 
      
    R = diag(diag(R));                  % R diagonal
    
    loglik = loglik_t;
    
    LL = [LL loglik];
    num_iter;
    
    initx = x1sum;
    initV = (P1sum - initx*initx');
    
    converged = em_converged(loglik, previous_loglik, thresh,1);
    
    previous_loglik = loglik;
    

end



[xitt,xittm,Ptt,Pttm,loglik_t]=K_filter(initx,initV,x,A,C,R,Q);
[xsmooth, Vsmooth, VVsmooth]=K_smoother(A,xitt,xittm,Ptt,Pttm,C,R);

chi = xsmooth'*C'*diag(Wx) + kron(ones(T,1),Mx);

F_hat =  xsmooth(1:r,:)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [beta, gamma, delta, gamma1, gamma2, x1, V1, loglik_t, xsmooth] = ...
    Estep(y, A, C, Q, R, initx, initV)

% This function computes the (expected) sufficient statistics for a single Kalman filter sequence.
%
% y is the observable and x the hidden state

% INPUTS
% y(:,t) - the observation at time t
% A - the system matrix
% C - the observation matrix 
% Q - the system covariance 
% R - the observation covariance
% initx - the initial state (column) vector 
% initV - the initial state covariance 

% OUTPUTS: the expected sufficient statistics, i.e.
% beta = sum_t=1^T (x_t * x'_t-1)
% gamma = sum_t=1^T (x_t * x'_t)  
% delta = sum_t=1^T (y_t * x'_t)
% gamma1 = sum_t=2^T (x_t-1 * x'_t-1)
% gamma2 = sum_t=2^T (x_t * x'_t)
% x1  expected value of the initial state
% V1  variance of the initial state
% loglik value of the loglikelihood
% xsmooth expected value of the state

[os T] = size(y);
ss = length(A);
 
                                                                            %   xsmooth = y;
                                                                            %   Vsmooth = zeros(ss, ss, T); % no uncertainty about the hidden states
                                                                            %   VVsmooth = zeros(ss, ss, T);
                                                                            %   loglik = 0;

% use the Kalman smoother to compute 
% xsmooth = E[X(:,t) | y(:,1:T)]
% Vsmooth = Cov[X(:,t) | y(:,1:T)]
% VVsmooth = Cov[X(:,t), X(:,t-1) | y(:,1:T)] t >= 2
% loglik = sum{t=1}^T log P(y(:,t))
[xitt,xittm,Ptt,Pttm,loglik_t]=K_filter(initx,initV,y',A,C,R,Q);
[xsmooth, Vsmooth, VVsmooth]=K_smoother(A,xitt,xittm,Ptt,Pttm,C,R);

% compute the expected sufficient statistics
delta = zeros(os, ss);
gamma = zeros(ss, ss);
beta = zeros(ss, ss);
for t=1:T
    delta = delta + y(:,t)*xsmooth(:,t)';
    gamma = gamma + xsmooth(:,t)*xsmooth(:,t)' + Vsmooth(:,:,t);
    if t>1 beta = beta + xsmooth(:,t)*xsmooth(:,t-1)' + VVsmooth(:,:,t); end
end
gamma1 = gamma - xsmooth(:,T)*xsmooth(:,T)' - Vsmooth(:,:,T);
gamma2 = gamma - xsmooth(:,1)*xsmooth(:,1)' - Vsmooth(:,:,1);

x1 = xsmooth(:,1);
V1 = Vsmooth(:,:,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xitt,xittm,Ptt,Pttm,loglik]=K_filter(initx,initV,x,A,C,R,Q)
% INPUTS
% x(:,t) - the observation at time t
% A - the system matrix
% C - the observation matrix 
% Q - the system covariance 
% R - the observation covariance
% initx - the initial state (column) vector 
% initV - the initial state covariance 
% OUTPUT:
% xittm = E[X(:,t) | y(:,1:t-1)]
% Pttm = Cov[X(:,t) | y(:,1:t-1)]
% xitt = E[X(:,t) | y(:,1:t)]
% Ptt = Cov[X(:,t) | y(:,1:t)]
%loglik - value of the loglikelihood

[T,N]=size(x);
r=size(A,1);

y=x';


xittm=[initx zeros(r,T)];
xitt=zeros(r,T);

Pttm=zeros(r,r,T);
Pttm(:,:,1)=initV;
Ptt=zeros(r,r,T);

for j=1:T
    
    L=inv(C*Pttm(:,:,j)*C'+R);
    nv=size(Pttm,1);
    xitt(:,j)=xittm(:,j)+Pttm(:,:,j)*C'*L*(y(:,j)-C*xittm(:,j));
    Ptt(:,:,j)=Pttm(:,:,j)-Pttm(:,:,j)*C'*L*C*Pttm(:,:,j); 
    xittm(:,j+1)=A*xitt(:,j);
    Pttm(:,:,j+1)=A*Ptt(:,:,j)*A'+Q;
    lik(j)=((2*pi)^(-N/2))*(abs((det(C*Pttm(:,:,j)*C'+R)))^(-.5))*...
        exp(-1/2*(y(:,j)-C*xittm(:,j))'*L*(-1/2*(y(:,j)-C*xittm(:,j))));
    
    
    e = y(:,j) - C*xittm(:,j); % error (innovation)
    n = length(e);
    ss = length(A);
    d = size(e,1);
    S = C*Pttm(:,:,j)*C' + R;
    GG = C'*diag(1./diag(R))*C;
    Sinv = inv(S);
    
    %%%%%%%%%%%%%%%%%%%%%%
    
    detS = prod(diag(R))*det(eye(ss)+Pttm(:,:,j)*GG);
    denom = (2*pi)^(d/2)*sqrt(abs(detS));
    mahal = sum(e'*Sinv*e,2);
    logl(j) = -0.5*mahal - log(denom);
    
end

loglik=sum(logl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xitT,PtT,PtTm]=K_smoother(A,xitt,xittm,Ptt,Pttm,C,R)
% INPUTS
% y(:,t) - the observation at time t
% A - the system matrix
% xittm = E[X(:,t) | y(:,1:t-1)]
% Pttm = Cov[X(:,t) | y(:,1:t-1)]
% xitt = E[X(:,t) | y(:,1:t)]
% Ptt = Cov[X(:,t) | y(:,1:t)]
% C - the observation matrix 
% R - the observation covariance

% OUTPUT:
% xitT = E[X(:,t) | y(:,1:T)]
% PtT = Cov[X(:,t) | y(:,1:T)]
% PtTm = Cov[X(:,t+1),X(:,t) | y(:,1:T)]

[T]=size(xitt,2);
r=size(A,1);
Pttm=Pttm(:,:,1:end-1);
xittm=xittm(:,1:end-1);
J=zeros(r,r,T);


for i=1:T-1
    J(:,:,i)=Ptt(:,:,i)*A'*inv(Pttm(:,:,i+1));
end

for i=1:T
    L(:,:,i)=inv(C*Pttm(:,:,i)*C'+R);
    K(:,:,i)=Pttm(:,:,i)*C'*L(:,:,i);
end


xitT=[zeros(r,T-1)  xitt(:,T)];
PtT=zeros(r,r,T);
PtTm=zeros(r,r,T);
PtT(:,:,T)=Ptt(:,:,T);
PtTm(:,:,T)=(eye(r)-K(:,:,T)*C)*A*Ptt(:,:,T-1);

for j =1:T-1
    
    xitT(:,T-j)= xitt(:,T-j)+J(:,:,T-j)*(xitT(:,T+1-j)-xittm(:,T+1-j));
    PtT(:,:,T-j)=Ptt(:,:,T-j)+J(:,:,T-j)*(PtT(:,:,T+1-j)-Pttm(:,:,T+1-j))*J(:,:,T-j)';
    
    
end


for j =1:T-2
    PtTm(:,:,T-j)=Ptt(:,:,T-j)*J(:,:,T-j-1)'+J(:,:,T-j)*(PtTm(:,:,T-j+1)-A*Ptt(:,:,T-j))*J(:,:,T-j-1)';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XC = center(X)
%CENTER XC = center(X)
%	Centers each column of X.

%	J. Rodrigues 26/IV/97, jrodrig@ulb.ac.be

[T n] = size(X);
XC = X - ones(T,1)*(sum(X)/T); % Much faster than MEAN with a FOR loop
