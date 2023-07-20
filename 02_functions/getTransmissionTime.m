function [time_tx,dtS] = getTransmissionTime(time_tx_RAW,Eph)
% ======================================================================= %
% Function to compute the transmission time of the satellite signal using
% clock corrections from broadcast ephemerides
% (based on 'GPSW Interface Specification IS-GPS-200 Revision E', p.86)
% ----------------------------------------------------------------------- %
% Input:
%   time_tx_RAW ... approximate transmission time (SOW)        double[1x1]
%           Eph ... ephemerides of satellite                   double[35x1]
%
% Output:
%       time_tx ... transmission time (SOW)                    double[1x1]
%           dtS ... clock correction (sec)                     double[1x1]
% ----------------------------------------------------------------------- %
% @author: Florian Zimmermann, Christian Eling   
% @date: 22.01.2014
% @mail: zimmermann@igg.bonn.de, eling@igg.uni-bonn.de
% ======================================================================= %

% ----- compute relativistic clock correction --------------------------- %
dtRel = getRelativisticClockCorr(time_tx_RAW,Eph);

% ----- compute clock correction for approximate transmission time ------ %
dtS = getSatClockCorr(time_tx_RAW,Eph);
dtS = dtS + dtRel;

% ----- recompute clock correction for corrected transmission time ------ %
dtS = getSatClockCorr(time_tx_RAW - dtS,Eph);
dtS = dtS + dtRel;

% ----- compute transmission time --------------------------------------- %
time_tx = time_tx_RAW - dtS;

return;



function dtRel = getRelativisticClockCorr(time,eph)
% ======================================================================= %
% Function to compute the relativistic correction term for satellite clocks
% using broadcast ephemerides
% (based on 'GPSW Interface Specification IS-GPS-200 Revision E', p.86)
% ----------------------------------------------------------------------- %
% Input:
%   time ... signal receiving time (SOW)             double[1x1]
%    eph ... ephemerides of satellite                double[35x1]
%
% Output:
%    dtS ... relativistic correction (sec)           double[1x1]
% ----------------------------------------------------------------------- %
% @author: Florian Zimmermann   
% @date: 22.01.2014
% @mail: zimmermann@igg.bonn.de
% ======================================================================= %

% geocentric gravitational constant
mu = 3.986005e14; % ref: Hoffmann-Wellenhoff

% ----- get parameters from ephemerides --------------------------------- %
delta_n = eph(13);              % mean motion difference from computed value
M_null = eph(14);               % mean anomaly at reference time
e = eph(16);                    % eccentricity
sqrt_a = eph(18);               % square root of semi-major axis
toe = eph(19);                  % reference time ephemeries
tgd = eph(33);                  % group delay 

% ----- compute relativistic correction --------------------------------- %
a = sqrt_a*sqrt_a;              % semi-major axis
n = sqrt(mu/(a^3))+delta_n;     % corrected mean motion
tk = check_sow(time-toe);       % time from ephemeries reference epoch

Mk = M_null+n*tk;               % mean anomaly

try
    Mk = rem(Mk+2*pi,2*pi);
catch
    disp 'ERROR Clock'
    keyboard;
end

E_null = Mk;                    % eccentric anomaly (solved by iteration)
differenz = 1;
while differenz>10^-12
    Ek = Mk+e*sin(E_null);
    diff = rem(Ek-E_null,2*pi);
    differenz = abs(diff);
    E_null = Ek;
end
Ek = rem(Ek+2*pi,2*pi);   
    
dtRel = -4.442807633e-10 * e * sqrt_a * sin(Ek);
dtRel = dtRel - tgd;            % relativistic correction and group delay

return;


function dtS = getSatClockCorr(time,eph)
% ======================================================================= %
% Function to compute the satellite clock correction using broadcast
% ephemerides
% (based on 'GPSW Interface Specification IS-GPS-200 Revision E', p.86)
% ----------------------------------------------------------------------- %
% Input:
%   time ... signal receiving time (SOW)             double[1x1]
%    eph ... ephemerides of satellite                double[35x1]
%
% Output:
%    dtS ... clock correction (sec)                  double[1x1]
% ----------------------------------------------------------------------- %
% @author: Florian Zimmermann   
% @date: 22.01.2014
% @mail: zimmermann@igg.bonn.de
% ======================================================================= %

% ----- get parameters from ephemerides --------------------------------- %
% polynomial coefficients
af0 = eph(8); % SV clock bias
af1 = eph(9); % SV clock drift
af2 = eph(10); % SV clock drift rate
% reference epoch of clock
toc = eph(19); 

% ----- compute clock correction ---------------------------------------- %
dt = check_sow(time - toc);
dtS = af0 + (af1 + af2*dt)*dt;

return; 

