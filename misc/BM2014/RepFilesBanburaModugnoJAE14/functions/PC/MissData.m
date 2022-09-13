%______________________________________________________________________
function [y,Z,G,c1,L]  = MissData(y,S);
%______________________________________________________________________
% PROC missdata                                                        
% PURPOSE: eliminates the rows in y & matrices Z, G that correspond to     
%          missing data (NaN) in y                                                                                  
% INPUT    y             vector of observations at time t  (n x 1 )    
%          S             KF system matrices             (structure)
%                        must contain Z & G
% OUTPUT   y             vector of observations (reduced)   (# x 1)     
%          Z G           KF system matrices     (reduced)   (# x ?)     
%          L             To restore standard dimensions     (n x #)     
%                        where # is the nr of available data in y
%______________________________________________________________________
  ix = ~isnan(y);
  e  = eye(size(y,1));
  L  = e(:,ix);

  y  =    y(ix);
  Z  =  S.Z(ix,:);  
  G  =  S.G(ix,:);
  c1 = S.c1(ix,:);
