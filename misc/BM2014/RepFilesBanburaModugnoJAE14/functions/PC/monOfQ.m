function m = monOfQ(m)
% function m = monOfQ(m)
% Delivers the month within the quarter 
% INPUT  m                      integer
% OUTPUT m                   in {1,2,3}
%______________________________________

  m = mod(m,3);
  m = m + 3*(m==0);
