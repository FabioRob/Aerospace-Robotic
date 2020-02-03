function [ v_inf_dep , v_dep , v_inf_arr , v_arr , h_H , e_H , a_H , TOF_H  ] = Compute_Direct_Hohmann (planet_id_dep , planet_id_arr  )

% -----------------------------------------------------------------------%
%
% Compute_Direct_Hohmann function computes direct Hohmann tranfer elements from two planets
%
% Arguments :
%
% planet_id_dep - departure planet,from 1 to 7 (from Mercury to Uranus)
% planet_id_arr - arrival planet, from 1 to 7 (from Mercury to Uranus)
%
% Outputs :
%
% v_inf_dep     - hyperbolic excess speed of departure hyperbola
% v_dep         - departure velocity of the spacecraft
% v_inf_arr     - hyperbolic excess speed of arrival hyperbola
% v_arr         - arrival velocity of the spacecraft
% h_H           - specific angular momentum of Hohmann transfer ellipse(Km^2/s)
% e_H           - eccentricity of Hohmann transfer ellipse
% a_H           - semiajor axis of the elliptical Hohmann transfer(Km)
% TOF_H         - period of the Hohmann transfer(sec)
%
% -----------------------------------------------------------------------%

muSun= 132712439935; % Sun's gravitational parameter(km^3/s^2)

[J2000_coe_dep ,~ ] = planetary_elements(2,planet_id_dep) ;
[J2000_coe_arr ,~ ] = planetary_elements(2,planet_id_arr) ;

r_dep = J2000_coe_dep(1); % departure planet orbit's radius(Km)
r_arr = J2000_coe_arr(1); % arrival planet orbit's radius(Km)

% Specific angular momentum of Hohmann transfer ellipse(Km^2/s)
h_H = sqrt ( 2*muSun * r_dep * r_arr / (r_dep+r_arr));
% Eccentricity of Hohmann transfer ellipse
e_H = (r_arr - r_dep) / (r_arr + r_dep);

% Heliocentric velocity of departure planet
v_pl_dep = sqrt(muSun/r_dep); 

% Velocity of the spacecraft at Earth departure
v_dep = h_H / r_dep;

% Spacecraft's velocity relative to departure planet at exit from its SOI
% to put it in Hohmann heliocentric
% transfer to arrival planet (it's the hyperbolic excess speed
% of the departure hyperbola)

v_inf_dep = v_dep - v_pl_dep ;

% Heliocentric velocity of arrival planet
v_pl_arr = sqrt(muSun/r_arr); 

% Velocity of the spacecraft on arrival
v_arr = h_H / r_arr;

% Spacecraft's velocity relative to Uranus at entrance in Uranus' SOI
% to put in in orbit around the planet (it's the hyperbolic excess speed
% of the arrival hyperbola)

v_inf_arr = v_pl_arr - v_arr ;

% Semiajor axis of the elliptical Hohmann transfer
a_H = (r_dep+r_arr)/2;

% Period of the Hohmann transfer(sec)
T_H = (2*pi / sqrt(muSun) ) * a_H ^1.5 ;

% Time of flight (TOF) from departure from Earth to arrival to Uranus(sec)
TOF_H = T_H / 2;