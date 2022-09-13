% Program to simulate the data and estimate the parameters (and factors) 
% of the dynamic factor model as in Section 3 of 
% Marta Banbura and Michele Modugno 
% "MAXIMUM LIKELIHOOD ESTIMATION OF FACTOR MODELS ON DATA SETS WITH 
% ARBITRARY PATTERN OF MISSING DATA."

clc;
clear;
close all;

%% set the simulation parameters
alpha = .7;       %% ar on factors
a     = 0;       %% ar on idio
b     = 0;       %% cross correlation of idio shocks
r     =  1;       %% # of dynamic factors
s     =  0;       %% # of lags of the dynamic factors


%% set the estimation parameters
r_hat        = r*(s+1); %% # of static factors
q_hat        = r ;      %% # of dynamic factors
p_hat        = 1;       %% # length of ar filter on common factors

P.p = p_hat;
P.r = r_hat;
P.max_iter = 500;    %% max # of iterations for ML estimation

% simulate data from factor model
[X,FF,tmp,beta] = sim_DGR(100,50,alpha,a,b,r,s,0.5*ones(50,1));
XTrue = X;
%missing data (10% of sample)
idx = randsample(numel(X),ceil(numel(X)*.1));
X(idx(1:ceil(numel(X)*.1)))  = nan;

% estimate the model and factors
Res_ml = EM_DFM(X,P);
Fml = Res_ml.F;

% rescale the factor estimate
b = (Fml'*Fml)\(Fml'*FF);
Fmlresc = Fml*b;

% plot the true and the estimated factor
figure;plot([FF Fmlresc]); legend('True factor','Estimated factor')

% Input data with missing values replaced by estimates
Xest = Res_ml.y;
% plot the true and the estimated data (first column)
figure;plot(1:100,XTrue(:,1),'b',1:100, Xest(:,1),'r'); legend('True data','Estimated data')
