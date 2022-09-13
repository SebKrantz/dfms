clear
clc

% Set Base_Dir
Base_Dir  = cd;
addpath(genpath(Base_Dir));


% start of the estimation sample
P.StartEst = [1993 1];
% first forecast produced in 
P.StartEv  = [1999 10];
% last forecast produced in 
P.EndEv    = [2008 3];
% forecast horizon
P.fcstH_Q = 1;


P.DataFile = 'Data\Data2009-10-15';
P.Legend = 'Data\Legend';

% data set
P.Model = 'Medium';


% Best ex-post parameterisations
% with AR(1) idiosyncratic component
P.idio = '_idio';
P.p = 2;P.r = [5];   funPRT_ML(P);
% with serially uncorrelated idiosyncratic component
P.idio = '';
P.p = 2;P.r = [3];   funPRT_ML(P);


% Remaining parameterisations for forecast averages
P.idio = '_idio';
P.p = 1;
P.r = [1];   funPRT_ML(P);P.r = [2];   funPRT_ML(P);P.r = [3];   funPRT_ML(P);
P.r = [4];   funPRT_ML(P);P.r = [5];   funPRT_ML(P);
P.p = 2;
P.r = [1];   funPRT_ML(P);P.r = [2];   funPRT_ML(P);P.r = [3];   funPRT_ML(P);
P.r = [4];   funPRT_ML(P);

P.idio = '';
P.p = 1;
P.r = [1];   funPRT_ML(P);P.r = [2];   funPRT_ML(P);P.r = [3];   funPRT_ML(P);
P.r = [4];   funPRT_ML(P);P.r = [5];   funPRT_ML(P);
P.p = 2;
P.r = [1];   funPRT_ML(P);P.r = [2];   funPRT_ML(P);
P.r = [4];   funPRT_ML(P);P.r = [5];   funPRT_ML(P);

