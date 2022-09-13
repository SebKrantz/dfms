clear
clc

% Set Base_Dir
Base_Dir  = cd;
addpath(genpath(Base_Dir));


P.DataFile = 'Data\Data2009-10-15';
P.Legend = 'Data\Legend';

% forecast series
P.SerNews  = 'GDP';
% forecast reference period
P.Qnews    = [2008 12];
% first forecast done in
P.StartEv  = [2008 7];
% last forecast done in
P.EndEv    = [2009 1];
% start of the estimation sample
P.StartEst = [1993 1];

% with AR(1) idiosyncratic component 
P.idio = '_idio'; 

P.r = 4; P.p = 2;
P.Model = 'Small';funNews_ML(P);
P.r = 5; P.p = 2;
P.Model = 'Medium';funNews_ML(P);
P.r = 5; P.p = 2;
P.Model = 'Large';funNews_ML(P);


