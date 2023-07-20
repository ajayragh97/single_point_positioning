function satPos=GetEphemeries(sow,eph)
% ======================================================================= %
% Function to determine ECEF-coordinates from broadcast ephemeries
% ----------------------------------------------------------------------- %
% Input:
%
%   sow ... time (SOW)                                    double[1x1]
%   eph ... broadcast ephemeries for single satellite     double[35x1]
%
% Output:
%   satPos ... satellite position [m]                     double[3x1]
%      rel ... relativistic time correction               double[1x1]
%      ura ... user range accuracy (std ephemeries) [m]   double[1x1]
% ----------------------------------------------------------------------- %
% @author: Christian Eling
% @revised by: Florian Zimmermann
% @date: 05.12.2013
% @mail: eling@igg.uni-bonn.de, zimmermann@igg.bonn.de
% ======================================================================= %

% geocentric gravitational constant
mu = 3.986005e14; % ref: Hoffmann-Wellenhoff
% speed of light
c = 299792458;
% earth rotation rate [rad*s^-1]
rot_rate = 7.2921151467e-5;

% ----- get parameters from eph ----------------------------------------- %
% first row of navigation message
af0 = eph(8);             % SV clock bias
af1 = eph(9);             % SV clock drift
af2 = eph(10);            % SV clock drift rate

% second row of navigation message
Crs = eph(12);            % amplitude of the sine harmonic correction term to the orbit radius
delta_n = eph(13);        % mean motion differemce from computed value
M_null = eph(14);         % mean anomaly at reference epoch

% third row of navigation message
Cuc = eph(15);            % amplitude of the cosine harmonic correction term to the argument of latitude
e = eph(16);              % eccentricity
Cus = eph(17);            % amplitude of the sine harmonic correction term to the argument of latitude
sqrt_a = eph(18);         % square root of the semi-major axis

% fourth row of navigation message
toe = eph(19);            % reference time ephemeries
Cic = eph(20);            % amplitude of the cosine harmonic correction term to the angle of inclination
Omega = eph(21);          % longitude of ascending node of orbit plane at weekly epoch
Cis = eph(22);            % amplitude of the sine harmonic correction term to the angle of inclination

% fifth row of navigation message
i_null = eph(23);         % inclination angle at reference time
Crc = eph(24);            % amplitude of the cosine harmonic correction term to the orbit radius
omega = eph(25);          % argument of perigee
Omega_dot = eph(26);      % rate of right ascension

% sixth row of navigation message
i_dot = eph(27);          % rate of inclination angle

% seventh row of navigation message
sva = eph(31);            % sv accuracy
TGD = eph(33);            % group delay

% eight row of navigation message
toc = eph(35);            % reference time clock

% ----- determination of satellite coordinates -------------------------- %
% based on 'Hofmann-Wellenhof p.432' and 'rtklib Manual'

a = sqrt_a*sqrt_a;          % semi-major axis
n = sqrt(mu/(a^3))+delta_n; % corrected mean motion

tk = check_sow(sow-toe);    % time from ephemeries reference epoch

Mk = M_null+n*tk;           % mean anomaly         
Mk = rem(Mk+2*pi,2*pi);

E_null = Mk;                % eccentric anomaly (solved by iteration)
differenz = 1;
while differenz>10^-12
    Ek = Mk+e*sin(E_null);
    diff = rem(Ek-E_null,2*pi);
    differenz = abs(diff);
    E_null = Ek;
end
Ek = rem(Ek+2*pi,2*pi);


v = atan2(sqrt(1-e^2)*sin(Ek), cos(Ek)-e); % true anomaly

phi_k = v+omega;            % argument of latitude
phi_k = rem(phi_k,2*pi);

d_phi_k = Cus*sin(2*phi_k)+Cuc*cos(2*phi_k); % argument of latitude correction
d_r_k = Crs*sin(2*phi_k)+Crc*cos(2*phi_k);   % radius correction
d_i_k = Cis*sin(2*phi_k)+Cic*cos(2*phi_k);   % inclination correction

u_k = phi_k+d_phi_k;         % corrected argument of latitude
r_k = a*(1-e*cos(Ek))+d_r_k; % corrected radius
i_k = i_null+i_dot*tk+d_i_k; % corrected inclination

Omega_k = Omega+(Omega_dot-rot_rate)*tk-rot_rate*toe; % corrected longitude of ascending node
Omega_k = rem(Omega_k+2*pi,2*pi);

% position in orbital plane
x_p = r_k*cos(u_k);
y_p = r_k*sin(u_k);
% ECEF-coordinates
x_s = x_p*cos(Omega_k)-y_p*cos(i_k)*sin(Omega_k);
y_s = x_p*sin(Omega_k)+y_p*cos(i_k)*cos(Omega_k);
z_s = y_p*sin(i_k);
satPos = [x_s;y_s;z_s];

