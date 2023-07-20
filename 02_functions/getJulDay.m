function jd = igg_getJulDay(y,m,d,h)
% ======================================================================= %
% Function to convert a date to julian date
% ----------------------------------------------------------------------- %
% Input:
%   y ... year (e.g. 2013)           int[1x1]
%   m ... month                      int[1x1]
%   d ... day                        int[1x1]
%   h ... hour and fraction hereof   double[1x1]
%
% Output:
%   jd ... julian date               double[1x1]
%
% The conversion is only valid in the time span from
% March, 1, 1900 to February, 28, 2100
% ----------------------------------------------------------------------- %
% @author: C. Eling
% @revised by: Florian Zimmermann
% @date: 05.12.2013
% @mail: zimmermann@igg.bonn.de
% ======================================================================= %

if m <= 2, y = y-1; m = m+12; end

jd = floor(365.25*(y+4716))+floor(30.6001*(m+1))+d+h/24-1537.5;
