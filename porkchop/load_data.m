global mu ephname km ip1 ip2 revmax rtd phem_populated

% taken 'lnxp1900p2053.421' frome here: ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/de421/
% and renamed to 'de421.bin'
% more info here: ftp://ssd.jpl.nasa.gov/pub/eph/planets/README.txt
ephname = 'de421.bin';

ip1 = 3; % launch planet is Earth
ip2 = 2; % arrival planet is Venus

% restrict to type I and II trajectories
revmax = 0;

% gravitational constant of the sun (km^3/sec^2)
mu = 132712441933.0;

km = 1; % km/s instead of AU/day

% rads to degree
rtd = 180.0 / pi;

% if variable 'phem_populated' not exist or it's empty
if exist('phem_populated','var') == 0 || length(phem_populated) == 0
    phem_populated = -1;
end