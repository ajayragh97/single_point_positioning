function [X, clkCorr] = satellitePosition(time,ephemeris,satellite)

% ======================================================================= %
% Function to determine ECEF-coordinates from broadcast ephemeries
% ----------------------------------------------------------------------- %
% Input:
%
%   time ... Gps seconds of week                          double[1x1]
%   eph ... broadcast ephemeries for all visible 
%              satellites and times                       double[35xn]
%
% Output:
%        X ... satellite position [m]                     double[3x1]
%  clkCorr ... Satellitenuhrenfehler                      double[1x1]
% ----------------------------------------------------------------------- %
% @author: Christian Eling
% @date: 10.11.2014
% @mail: eling@igg.uni-bonn.de, zimmermann@igg.bonn.de
% ======================================================================= %

% satellite orbits in the navigation file are mostly available for several
% points in time. The nearst parameter set is used for the following orbit
% determination
icol = find_eph(ephemeris,satellite,time);

if icol<1 || isempty(icol)
    disp 'Error Ephemeris....'
end

% calculation of clock corrections for satellite clock 
[timeCorr, clkCorr] = getTransmissionTime(time,ephemeris(:,icol));

 % calculation of satellite position
X = GetEphemeries(timeCorr, ephemeris(:,icol));