function [lat, lon, height] = convertECEFtoGRS80(x,y,z)
% Convert Earth-centered Earth-fixed (ECEF) coordinates 
% to geodetic coordinates(GRS80).
a = 6378137; % GRS80
f = 1 / 298.257222101; % GRS80
b = a * (1 - f);
e_2 = (a.^2 - b.^2) / a.^2;
e_prime_power_2 = (a.^2 - b.^2) / b.^2;
p = sqrt(x.^2 + y.^2);
h = 0;
lat = atan2(z * (1 + e_prime_power_2), p);

difference_latitude = 1;
difference_height = 1;
while any(difference_latitude(:) > 0.00001) || any(difference_height(:) > 0.00001)
    starting_lat = lat;
    starting_h = h;
    n = a./ sqrt(1 - e_2 * sin(lat).^2);
    h = p./ cos(lat) - n;
    lat = atan2(z.* (n + h), p.* ((n./ (1 + e_prime_power_2)) + h));
    difference_latitude = abs(starting_lat - lat);
    difference_height = abs(starting_h - h);
end

lat = lat;
lon = atan2(y,x);
height = h;

end