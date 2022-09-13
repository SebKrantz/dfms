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


% Specification of the series list
P.Model = 'BaRu'; 

% Parameterisation
P.r = 5; P.q = 2 ;P.p = 2;

% Running the out-of-sample forecast
funPRT_PC(P);
