function [J2000_coe, rates] = planetary_elements(assumption,planet_id)
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
%
% This function extracts a planet’s J2000 orbital elements and
% centennial rates from Table 8.1.
%
% if assumption=1 it assumes real orbits
% if assumption=2 it assumes circular and coplanar orbits(for Direct Hohmann transfer)
% if assumption=3 it assumes coplanar orbits(for design with flyby)

% planet_id - 1 through 9, for Mercury through Pluto
%
% J2000_elements - 9 by 6 matrix of J2000 orbital elements for
% the nine planets Mercury through Pluto. The
% columns of each row are:
% a = semimajor axis (Km)
% e = eccentricity  -- zero if assumption=2
% i = inclination (degrees)  -- zero if assumption=2 or 3
% RA = right ascension of the ascending -- zero if assumption=2 or 3
% node (degrees)
% w_hat = longitude of perihelion (degrees) -- zero if assumption=2
% L = mean longitude (degrees)
%

% cent_rates - 9 by 6 matrix of the rates of change of the
% J2000_elements per Julian century (Cy).
% Using ''dot'' for time derivative, the
% columns of each row are:
% a_dot (AU/Cy)
% e_dot (1/Cy)
% i_dot (arcseconds/Cy)
% RA_dot (arcseconds/Cy)
% w_hat_dot (arcseconds/Cy)
% Ldot (arcseconds/Cy)
%
% J2000_coe - row vector of J2000_elements corresponding
% to ''planet_id'', with au converted to km
% rates - row vector of cent_rates corresponding
% to ''planet_id'', with au converted to km
% and arcseconds converted to degrees
%
% au - astronomical unit (km)
%
% User M-functions required: none
% ------------------------------------------------------------
if assumption==1 % real orbits
    J2000_elements = ...
    [ 0.38709893 0.20563069 7.00487 48.33167 77.45645 252.25084
    0.72333199 0.00677323 3.39471 76.68069 131.53298 181.97973
    1.00000011 0.01671022 0.00005 -11.26064 102.94719 100.46435
    1.52366231 0.09341233 1.85061 49.57854 336.04084 355.45332
    5.20336301 0.04839266 1.30530 100.55615 14.75385 34.40438
    9.53707032 0.05415060 2.48446 113.71504 92.43194 49.94432
    19.19126393 0.04716771 0.76986 74.22988 170.96424 313.23218
    30.06896348 0.00858587 1.76917 131.72169 44.97135 304.88003
    39.48168677 0.24880766 17.14175 110.30347 224.06676 238.92881];
    cent_rates = ...
    [ 0.00000066 0.00002527 -23.51 -446.30 573.57 538101628.29
    0.00000092 -0.00004938 -2.86 -996.89 -108.80 210664136.06
    -0.00000005 -0.00003804 -46.94 -18228.25 1198.28 129597740.63
    -0.00007221 0.00011902 -25.47 -1020.19 1560.78 68905103.78
    0.00060737 -0.00012880 -4.15 1217.17 839.93 10925078.35
    -0.00301530 -0.00036762 6.11 -1591.05 -1948.89 4401052.95
    0.00152025 -0.00019150 -2.09 -1681.4 1312.56 1542547.79
    -0.00125196 0.00002514 -3.64 -151.25 -844.43 786449.21
    -0.00076912 0.00006465 11.07 -37.33 -132.25 522747.90];
elseif assumption==2 % circular and coplanar orbits
    J2000_elements = ...
    [ 0.38709893 0.0 0.0 0.0 0.0 252.25084
    0.72333199 0.0 0.0 0.0 0.0 181.97973
    1.00000011 0.0 0.0 0.0 0.0 100.46435
    1.52366231 0.0 0.0 0.0 0.0 355.45332
    5.20336301 0.0 0.0 0.0 0.0 34.40438
    9.53707032 0.0 0.0 0.0 0.0 49.94432
    19.19126393 0.0 0.0 0.0 0.0 313.23218
    30.06896348 0.0 0.0 0.0 0.0 304.88003
    39.48168677 0.0 0.0 0.0 0.0 238.92881];
    cent_rates = ...
    [ 0.00000066 0.0 0.0 0.0 0.0 538101628.29
    0.00000092 0.0 0.0 0.0 0.0 210664136.06
    -0.00000005 0.0 0.0 0.0 0.0 129597740.63
    -0.00007221 0.0 0.0 0.0 0.0 68905103.78
    0.00060737 0.0 0.0 0.0 0.0 10925078.35
    -0.00301530 0.0 0.0 0.0 0.0 4401052.95
    0.00152025 0.0 0.0 0.0 0.0 1542547.79
    -0.00125196 0.0 0.0 0.0 0.0 786449.21
    -0.00076912 0.0 0.0 0.0 0.0 522747.90];
elseif assumption==3 % coplanar orbits
    J2000_elements = ...
    [ 0.38709893 0.20563069 0.0 0.0 77.45645 252.25084
    0.72333199 0.00677323 0.0 0.0 131.53298 181.97973
    1.00000011 0.01671022 0.0 0.0 102.94719 100.46435
    1.52366231 0.09341233 0.0 0.0 336.04084 355.45332
    5.20336301 0.04839266 0.0 0.0 14.75385 34.40438
    9.53707032 0.05415060 0.0 0.0 92.43194 49.94432
    19.19126393 0.04716771 0.0 0.0 170.96424 313.23218
    30.06896348 0.00858587 0.0 0.0 44.97135 304.88003
    39.48168677 0.24880766 0.0 0.0 224.06676 238.92881];
    cent_rates = ...
    [ 0.00000066 0.00002527 0.0 0.0 573.57 538101628.29
    0.00000092 -0.00004938 0.0 0.0 -108.80 210664136.06
    -0.00000005 -0.00003804 0.0 0.0 1198.28 129597740.63
    -0.00007221 0.00011902 0.0 0.0 1560.78 68905103.78
    0.00060737 -0.00012880 0.0 0.0 839.93 10925078.35
    -0.00301530 -0.00036762 0.0 0.0 -1948.89 4401052.95
    0.00152025 -0.00019150 0.0 0.0 1312.56 1542547.79
    -0.00125196 0.00002514 0.0 0.0 -844.43 786449.21
    -0.00076912 0.00006465 0.0 0.0 -132.25 522747.90];

    
else
    disp('first argument must be 1 for real ,2 for coplanar and circular, 3 for coplanar orbits')
    return
end
    
J2000_coe = J2000_elements(planet_id,:);
rates = cent_rates(planet_id,:);
%...Convert from AU to km:
au = 149597871;
J2000_coe(1) = J2000_coe(1)*au;
rates(1) = rates(1)*au;
%...Convert from arcseconds to fractions of a degree:
rates(3:6) = rates(3:6)/3600;

return