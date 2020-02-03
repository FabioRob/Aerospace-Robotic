function [coe, r, v, jd] = planet_elements_and_sv(assumption,planet_id, year, month, day, hour, minute, second)
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
%
% This function calculates the orbital elements and the state
% vector of a planet from the date (year, month, day)
% and universal time (hour, minute, second).
% if assumption=1 it assumes real orbits
% if assumption=2 it assumes circular and coplanar orbits(for Direct Hohmann transfer)
% if assumption=3 it assumes coplanar orbits(for design with flyby)

% mu - gravitational parameter of the sun (kmˆ3/sˆ2)
% deg - conversion factor between degrees and radians
% pi - 3.1415926...
%
% coe - vector of heliocentric orbital elements
% [h e RA incl w TA a w_hat L M E],
% where
% h = angular momentum (kmˆ2/s)
% e = eccentricity -- zero if assumption=2
% RA = right ascension (deg)  -- zero if assumption=2 or 3
% incl = inclination (deg)  -- zero if assumption=2 or 3
% w = argument of perihelion (deg) -- zero if assumption=2
% TA = true anomaly (deg)
% a = semimajor axis (km)
% w_hat = longitude of perihelion
% ( = RA + w) (deg)
% L = mean longitude ( = w_hat + M) (deg)
% M = mean anomaly (deg)
% E = eccentric anomaly (deg)
%
% planet_id - planet identifier:
% 1 = Mercury
% 2 = Venus
% 3 = Earth
% 4 = Mars
% 5 = Jupiter
% 6 = Saturn
% 7 = Uranus
% 8 = Neptune
% 9 = Pluto
%
% year - range: 1901 - 2099
% month - range: 1 - 12
% day - range: 1 - 31
% hour - range: 0 - 23
% minute - range: 0 - 60
% second - range: 0 - 60
%
% j0 - Julian day number of the date at 0 hr UT
% ut - universal time in fractions of a day
% jd - julian day number of the date and time
%
% J2000_coe - row vector of J2000 orbital elements from
% Table 8.1
% rates - row vector of Julian centennial rates from
% Table 8.1
% t0 - Julian centuries between J2000 and jd
% elements - orbital elements at jd
%
% r - heliocentric position vector
% v - heliocentric velocity vector
%
% User M-functions required: J0, kepler_E, sv_from_coe
% User subfunctions required: planetary_elements, zero_to_360
% ------------------------------------------------------------
mu=132712439935; % sun's gravitational parameter

deg = pi/180;
%...Equation 5.48:
j0 = J0(year, month, day);
ut = (hour + minute/60 + second/3600)/24;
%...Equation 5.47
jd = j0 + ut;
%...Obtain the data for the selected planet from Table 8.1:
[J2000_coe, rates] = planetary_elements(assumption,planet_id);
%...Equation 8.104a:
t0 = (jd - 2451545)/36525;

%...Equation 8.104b:
elements = J2000_coe + rates*t0;
a = elements(1);
e = elements(2);
%...Equation 2.61:
h = sqrt(mu*a*(1 - e^2));
%...Reduce the angular elements to within the range 0 - 360 degrees:
incl = elements(3);
RA = zero_to_360(elements(4));
w_hat = zero_to_360(elements(5));
L = zero_to_360(elements(6));
w = zero_to_360(w_hat - RA);
M = zero_to_360((L - w_hat));
%...Algorithm 3.1 (for which M must be in radians)
E = kepler_E(e, M*deg); % E in radians
%...Equation 3.10 (converting the result to degrees):
TA = zero_to_360...
(2*atan(sqrt((1 + e)/(1 - e))*tan(E/2))/deg);
coe = [h e RA incl w TA a w_hat L M E/deg]; % all angles in degrees
%coe = [h e RA*deg incl*deg w*deg TA*deg a w_hat*deg L*deg M*deg E]; % all angles in radians
%...Algorithm 4.2 (for which all angles must be in radians):
[r, v] = sv_from_coe([h e RA*deg incl*deg w*deg TA*deg],mu);
return
