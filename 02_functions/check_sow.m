function tt = check_sow(t)
% ======================================================================= %
% Function to repair over- and underflow of GPS time
% ----------------------------------------------------------------------- %
% Input:
%    t ... time (SOW)    double[1x1]
%
% Output:
%   tt ... checked time (SOW)   double[1x1]
% ----------------------------------------------------------------------- %
% @author: Christian Eling   
% @revised by: Florian Zimmermann
% @date: 05.12.2013
% @mail: eling@igg.uni-bonn.de, zimmermann@igg.bonn.de
% ======================================================================= %

% SOW for half week
half_week = 302400;
tt = t;

% check overflow
if t >  half_week
     tt = t-2*half_week; 
end

% check underflow
if t < -half_week
     tt = t+2*half_week;
end

return;