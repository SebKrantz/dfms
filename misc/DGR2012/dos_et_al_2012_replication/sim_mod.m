function [X,F,Lambda,R] = sim_mod(T,N,alpha,a,b,r,s,Lambda,R);
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
                                % set 0, for heterogenous signal ro noise ratio
                               
                                
if nargin < 8                 % if parameters of the model are not specified,
    Lambda = randn(N,r*(s+1)); % simulate the factor loading from a normal distribution
    if spherical               
        R = .5*ones(N,1); %% Relative variance of the idiosyncratic equal to 1/2 for all the cross-section
    else
        R = (rand(N,1)*.9+.1); %% Relative variance of the idiosyncratic equal to 1/2 for all the cross-section
    end;
end;

% simulate T+200 shocks from N+r normal distributions
temp = randn(T+200,N+r);
% construct the common shocks
u = temp(:,1:r);
% construct the idiosyncratic shocks with cross-correlation b
v = temp(:,r+1:end)*chol(toeplitz(b.^(0:1:N-1)));

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


C = var(F(end-T+1:end,:)*Lambda')'.*R./(1-R);

% Generate the idiosyncrati componet
e = filter(1,[1 -a],v/sqrt(1/(1-a^2)))*diag(sqrt(C));

X = F*Lambda'+e;

X = X(end-T+1:end,:);           % the simulated observations

F = F(end-T+1:end,:);           % the simulated factors

