function [X,F,X_f,R] = sim_DGR(T,N,alpha,a,b,r,s,R,Lambda);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Quasi Maximum Likelihood Approach for Large Approximate Dynamic Factor Models (2006a). CEPR DP n. 5724
%
% A two-step estimator for large approximate dynamic factor models based on Kalman filtering (2006b). Manuscript, ECARES-ULB.
%
% Catherine Doz, Universite' Cergy-Pontoise
% Domenico Giannone, Universite' Libre de Bruxelles and ECARES,
% Lucrezia Reichlin, European Central Bank, ECARES and CEPR
%
% Programs and manuscript available at: http://homepages.ulb.ac.be/~dgiannon/
%
% The programs have been organized and commented by Laura Coroneo and Michele Modugno from ULB-ECARES.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function simulates data from the model
%
% X = F*Lambda' + e
%
% 1) X is a T X N matrix of simulated data;
%
% 2) F T X r(s+1) matrix of static factors (dynamic factors and thier lags), 
%    r is the number of dynamic factors, 
%    s is the number of lagged dynamic factors 
%
% 3) Lambda N X r(s+1) matrix of the factor laodings 
% 
% 4) e T X N matrix of the idiosyncratic component
%
% INPUTS:
% T      - number of observations
% N      - number of cross sections
% alpha  - ar parameter on factors
% a      - ar parameter on the idiosyncratic components
% b      - cross correlation of idiosyncratic shocks
% r      - # of dynamic factors
% s      - # of lags of dynamic factors
% Lambda - matrix of loadings (this input is optional, 
%          if it is not given the function generates
%          a matrix of loadings)
% R      - Matrix of commonalities



spherical = 0;                  % set 1 to have a noise to signal ratio equal to 1 for all elements, 
                                
                               
                                
if nargin < 8 
    if spherical               
        R = .5*ones(N,1); %% Relative variance of the idiosyncratic equal to 1/2 for all the cross-section
    else
        R = (rand(N,1)*.8+.1); %% Relative variance of the idiosyncratic from U[0.1 0.9]
    end;
end


if nargin < 9 
    % if parameters of the model are not specified,
    Lambda = randn(N,r*(s+1)); % simulate the factor loading from a normal distribution
end;

% simulate T+200 shocks from N+r normal distributions
temp = randn(T+200,N+r);
% construct the common shocks
u = temp(:,1:r);
% construct the idiosyncratic shocks with cross-correlation b

v = temp(:,r+1:end);
if  b ~= 0
    v = v*chol(toeplitz(b.^(0:1:N-1)));
end

F = [];

% construct F the static factors: 
% from r dynamic factors with AR coefficient alpha and their lags.

for i = 1:r
    clear Ftemp
    for j = 1:s+1
        filt = zeros(s+1,1); filt(j)=1;
        Ftemp(:,j) = filter(filt,[1 -alpha],u(:,i)/sqrt(1/(1-alpha^2)));
    end;
    F = [F Ftemp];
end;
F = F-repmat(mean(F(end-T+1:end,:)),[T+200 1]);

C = var(F(end-T+1:end,:)*Lambda')'.*R./(1-R);

% Generate the idiosyncrati componet
e = filter(1,[1 -a],v/sqrt(1/(1-a^2)))*diag(sqrt(C));
% a_ols = diag(inv(diag(diag(e(1:end-1,:)'*e(1:end-1,:))))*diag(diag(e(1:end-1,:)'*e(2:end,:))));

X = F*Lambda'+e;

X = X(end-T+1:end,:);           % the simulated observations

F = F(end-T+1:end,:);           % the simulated factors

X_f = F*Lambda';
