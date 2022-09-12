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
% This is the main file that performs a Montecarlo evaluation of the QML and two steps estimates of the common factors.
% The output of this program are the measures of performance reported in the paper.
%                       
%
% It uses the following functions.
% sim_mod:           generates time series from the simulation model.
% DynFA:             extracts the unobservable factors using three different methods 
%                       QML:      Max Likelihood estimates using the Expectation Maximization (EM) algorithm 
%                                 (Doz, Giannone and Reichlin, 2006a) 
%                       TWO STEP: Principal components+Kalman filtering 
%                                 (Doz, Giannone and Reichlin, 2006b)
%                       PC:       principal components 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('----------------------------------------------------------------------------------------------------------------------------')

disp('replication files:')
disp(' ')
disp('A Quasi Maximum Likelihood Approach for Large Approximate Dynamic Factor Models (2012)')
disp('The Review of Economics and Statistics, MIT Press, vol. 94(4), pages 1014-1024, November.')
disp(' ')
disp('Catherine Doz, Universite Cergy-Pontoise')
disp('Domenico Giannone, Universite Libre de Bruxelles, ECARES and CEPR')
disp('Lucrezia Reichlin, London Business School and CEPR')
disp(' ')
disp('Programs are also available at: http://homepages.ulb.ac.be/~dgiannon/')
disp(' ')
disp('----------------------------------------------------------------------------------------------------------------------------')

clear;

REP = 30;                   %% Number of draws for the parameters

repetitions = REP^2;        %% Total number of montecarlo repetitions

NN = [10 25];            %% cross-sectional dimension  NN = [10 25 50 100] in the paper

TT = [50];              %% sample size                 TT = [50 100] in the paper



% sets the simulation parameters
alpha = .9;       %% ar on factors
a     = .5;       %% ar on idio
b     = .5;       %% cross correlation of idio shocks
r     =  3;       %% # of dynamic factors
s     =  0;       %% # of lags of the dynamic factors



% sets the estimation parameters
r_hat        = r*(s+1); %% # of static factors
q_hat        = r ;      %% # of dynamic factors
p_hat        = 1;       %% # length of ar filter on common factors
max_iter     = 2000;    %% max # of iterations for ML estimation


for j = 1:repetitions
    
    if mod(j,REP)==1;; %% generate the parameters every REP replications
        [XX,FF,Lambda,R] = sim_mod(TT(end),NN(end),alpha,a,b,r,s); % simulate the data and the parameters at first repetition
    else
        [XX,FF] = sim_mod(TT(end),NN(end),alpha,a,b,r,s,Lambda,R); % simulate the data given the parameters of the first repetition
    end;
    
    [XX,FF,Lambda] = sim_mod(TT(end),NN(end),alpha,a,b,r,s); 
    
    for jn = 1:length(NN)
        for jt = 1:length(TT)
            F = FF(1:TT(jt),:);
            X = XX(1:TT(jt),1:NN(jn));
            tic
            % estimates the common factors with maximum likelihood(F_hat), the
            % pricipal components (F_pca) and the two steps estimates (F_kal)
            [F_hat,F_pca,F_kal,num_iter(jt,jn,j)] = DynFA(X,q_hat,r_hat,p_hat,max_iter); 
            
            elapsed_time(jt,jn,j) = toc;
            
            % compute the measure of peformance of the diferent estimators
            PF_pc = F_pca*inv(F_pca'*F_pca)*F_pca';
            PF_fa = F_hat*inv(F_hat'*F_hat)*F_hat';
            PF_kf = F_kal*inv(F_kal'*F_kal)*F_kal';
            
            tr_fa(jt,jn,j) = trace(F'*PF_fa*F)/trace(F'*F);
            tr_pc(jt,jn,j) = trace(F'*PF_pc*F)/trace(F'*F);
            tr_kf(jt,jn,j) = trace(F'*PF_kf*F)/trace(F'*F);
            
        end;
    end;
    
    if mod(j,5)==0
        disp('-------------------------------------------')
        disp('-------------------------------------------')
        disp(['now running the ',num2str(j),'th repetition'])
        disp('Computational time in seconds')
        disp([NaN NN; 
              TT' mean(elapsed_time,3)])
        disp('iterations')
        disp([NaN NN; 
              TT' ceil(mean(num_iter,3))])
        disp('---------------------------------------------')
        disp('Trace statistics')
        disp('Max. likelihood');
        disp([NaN NN; 
              TT' mean(tr_fa,3)])
        disp('------------------------------');
        disp('FA/PC');
        disp([NaN NN; 
              TT' mean(tr_fa./tr_pc,3)])
        disp('------------------------------');
        disp('FA/KF');
        disp([NaN NN; 
              TT' mean(tr_fa./tr_kf,3)])
        disp('---------------------------------------------')
        disp('Rows indicate the sample size; Columns indicate the cross-sectional size')

    end
    
end;
