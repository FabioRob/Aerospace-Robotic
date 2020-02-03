function [mu,radius,SOI] = physical_data(planet_id)

% -----------------------------------------------------------------------%
%
% physical_data function returns gravitational parameter,radius and SOI
% radius of a planet
%
% Arguments :
%
% planet_id  -  planet,from 1 to 7 (from Mercury to Uranus)
% 
% Outputs :
%
% mu         - gravitational parameter of central body(Km^3/s^2)
% radius     - planet's radius
% SOI        - planet SOI's radius
%
% -----------------------------------------------------------------------%

physical = ...
[22030,324900,398600,42828,126686000,37931000,5794000; % gravitational parameters
2440,6052,6378,3396,71490,60270,25560; % planetary radii
112000,616000,925000,577000,48200000,54800000,51800000 ]; % SOI

mu = physical(1,planet_id) ;
radius = physical(2,planet_id);
SOI = physical(3,planet_id);