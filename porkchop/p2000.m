function [r, v] = p2000(ncent, ntarg, jdate)

% lunar/planetary state vector

% EME2000 coordinate system (jpl ephemeris)

% input

%  jdate = julian date
%  ntarg = "target" body

% output

%  r = position vector (kilometers)
%  v = velocity vector (kilometers/second)

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% call jpl ephemeris

result = jplephem (jdate, ntarg, ncent);

% load position and velocity vectors

r = result(1: 3);
v = result(4: 6);