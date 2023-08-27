function [E, N, h] = ell2utm32(lat, lon, height)
% convert ellipsoidal lat long h to utm32 
% parameter_grs80

G_0 = 111132.952547; %  m/degree
G_2 = -16038.5088; % m
G_4 = 16.8326; % m
G_6 = -0.0220; % m
a = 6378137.000; % m
b = 6356752.314; % m

L = rad2deg(lon);
B_degree = rad2deg(lat);
B = lat;

L_0 = 9;
diff_L = L - L_0;
E_0 = ((L_0 + 3)/ 6 + 30.5) * 10.^6;
m = 0.9996;

RHO = 180 / pi;
MU_p2 = ((a^2 - b^2) / b^2) * cos(B).^2;
C = (a^2) / b;
t = tan(B);
N_BAR = C / sqrt(1 + MU_p2);

G = G_0 * B_degree + G_2 * sin(2*B) + G_4 * sin(4*B) + G_6 * sin(6*B);
eq_1 = m / RHO * N_BAR * cos(B);
eq_3 = m / (6 * RHO.^3).* N_BAR * cos(B).^3 * (1 - t.^2 + MU_p2);
eq_5 = m / (120 * RHO.^5).* N_BAR * cos(B).^5 * (5 - 18*t.^2 + t.^4 + MU_p2.* (14 - 58.*(t.^2)));
eq_2 = m / (2 * RHO.^2).* N_BAR * cos(B).^2 * t;
eq_4 = m / (24 * RHO.^4).* N_BAR * cos(B).^4 * t.* (5 - t.^2 + 9 * MU_p2);
eq_6 = m / (720 * RHO.^6).* N_BAR * cos(B).^6 * t.* (61 - 58*t.^2 + t.^4);

%     E = E_0   + eq_1 * diff_L    + eq_3 * diff_L**3 + eq_5 * diff_L**5
E = E_0   + eq_1.* diff_L    + eq_3.* diff_L.^3 + eq_5.* diff_L.^5 - 32000000;
N = m.* G + eq_2.* diff_L.^2 + eq_4.* diff_L.^4 + eq_6.* diff_L.^6;

geoid_undulation = 47.020; % http://gibs.bkg.bund.de/geoid/gscomp.php?p=g
h = height - geoid_undulation;

end