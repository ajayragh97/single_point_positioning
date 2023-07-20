function sow = igg_getSOW(Date)
% ======================================================================= %
% Function to determine the "Second of Week" from a given Date
% ----------------------------------------------------------------------- %
% Input:
%   Date ... time (year,month,day,hour,min,sec)    double[6x1]
%
% Output:
%    sow ... converted time (SOW)                  double[1x1]
% ----------------------------------------------------------------------- %
% @author: C. Eling   
% @revised by: Florian Zimmermann
% @date: 05.12.2013
% @mail: zimmermann@igg.bonn.de
% ======================================================================= %

% compute julian date
jd = getJulDay(Date(1),Date(2),Date(3),Date(4)+Date(5)/60+Date(6)/3600);

% compute Seconds of Week from julian date
[~,sow] = getGpsTime(jd);
sow = round(sow*100)/100;
% sow = igg_checkSOW(sow);
