function [week,sec_of_week] = igg_getGpsTime(julday)
% ======================================================================= %
% Function to convert Julian Day number to GPS week and Seconds of Week 
% Seconds of Week are reckoned form Saturday midnight
% ----------------------------------------------------------------------- %
% Input:
%    julday ... julian date                    double[1x1]
%
% Output:
%          week ... GPS week                   double[1x1]
%   sec_of_week ... GPS seconds of week        double[1x1]
% ----------------------------------------------------------------------- %
% @author: C. Eling   
% @revised by: Florian Zimmermann
% @date: 05.12.2013
% @mail: zimmermann@igg.bonn.de
% ======================================================================= %

a = floor(julday+.5);
b = a+1537;
c = floor((b-122.1)/365.25);
e = floor(365.25*c);
f = floor((b-e)/30.6001);
d = b-e-floor(30.6001*f)+rem(julday+.5,1);
day_of_week = rem(floor(julday+.5),7);
week = floor((julday-2444244.5)/7);
% We add +1 as the GPS week starts at Saturday midnight
sec_of_week = (rem(d,1)+day_of_week+1)*86400;
