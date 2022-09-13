function [Xc, Jc] = outliers(X,type,options)
%__________________________________________________________________________
% function z = outliers(x,type,options)
% Identifies missing values (NaN) and outliers in series x
% Replaces NaN and outliers according to input TYPE 
%
% INPUT   X       Data matrix                 (nobs x k)
%        
%         type  option for correction         (integer)
%               0 -> replace outliers         with    NaN
%               1 -> replace outliers         with median 
%               2 -> replace outliers and NaN with median
%               3 -> replace outliers         with largest admiss value
%                            NaN              with median
%         options
%         k    [3]   Length of MA filter L^(-k) + .. 1 + .. L^k
%         c    [4]   Cutoff for outliers is c x interquintile distance
%                    of available data
%         r    [4]   Largest admissable value (r x interquintile dist)
%
% Replacement works in two steps: a) replace with median of the series
%                                 b) apply centred MA to replaced obs
%
% OUTPUT  Xc          Corrected series              (nobs x k)
%         Jc          Index of corrections          (nobs x k)
%                     =  0 if original data
%                     =  1 if outlier replaced
%                     = -1 if missing replaced
%_______________________________________________________________________________
% 
  if nargin<3 options=[]; end
  if ~isfield(options,   'k');  options.k  = 3        ; end
  if ~isfield(options,   'c');  options.c  = 4        ; end
  if ~isfield(options,   'r');  options.r  = options.c; end

  [n,m] = size(X);
  Xc   = nan*zeros(n,m);
  Jc   =     zeros(n,m);

%_______________________________________________________________________________
% Go!  
  for i = 1:m  
      x = X(:,i);
    % Make indices of NaN
      Jmiss = isnan(x);                              % index of all NaN
      
      if ~all(Jmiss)
      % Calculate median 
        xsort  = sort(x(Jmiss == 0));
        ns     = size(xsort,1);
        median = xsort(round(ns/2));
        xsign  = sign(x - median);
      
      % Make indices of outliers
      % Define outliers as > options.c x interquintile distance of data
        IQD    = abs(xsort(round(ns*1/5))-median);
        cutoff = options.c * IQD;
        Joutr  = abs(x-median) > cutoff;
      
      % Adjust  
        if type == 0                           % Replace outliers with NaN
           x(Joutr)  = NaN;

        elseif type == 1;                      % Replace outliers  
           x(Joutr)  = median;
           x_MA      = MAcentered(x,options.k);
           x(Joutr)  = x_MA(Joutr);

        elseif type == 2;                       % Replace NaN & outliers 
           x(Jmiss)  = median;          
           x(Joutr)  = median; 
           x_MA      = MAcentered(x,options.k);
           x(Joutr)  = x_MA(Joutr);
           x(Jmiss)  = x_MA(Jmiss);
      
        elseif type == 3                        % Replace outliers with LAV
           x(Jmiss)  = median;          
           x(Joutr)  = options.r * IQD * xsign(Joutr); 
           x_MA      = MAcentered(x,options.k);
           x(Joutr)  = x_MA(Joutr);
           x(Jmiss)  = x_MA(Jmiss);
           
        elseif type == 4
           Joutr(end-options.lag:end,:) = logical(0);
           x(Joutr)  = NaN;
           
        elseif type == 5;                       % Replace NaN
           x(Jmiss)  = median;          
           x_MA      = MAcentered(x,options.k);
           x(Jmiss)  = x_MA(Jmiss);
      
        end
        Jc(Joutr,i) =  1;
      end
      
      Xc(:,i)     =  x;
      Jc(Jmiss,i) = -1;

  end;
  
  
 
%_______________________________________________________________________________
function y = MAcentered(x,k);
%_______________________________________________________________________________
% Computes centered moving average y of order 2*k+1 of series x along columns
% x is expanded on both sides with 1st and last value respectively

  y  = zeros(size(x));
  x  = [x(1,:).*ones(k,1); x; x(end,:).*ones(k,1)];

  for i = k+1:size(x,1)-k
      y(i-k,:) = mean(x(i-k:i+k,:));
  end
 