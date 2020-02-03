function [] = plotPlanetOrbit (assumption , planet_id , date , col , linewidth)

% -----------------------------------------------------------------------%
%
% This function plots planet orbit real/coplanar and circular/coplanar(for
% assumption =1,2,3 respectevely )
%
% Arguments :
%
% assumption - 1,2,3 for,repectively, real,coplanar and circular, coplanar orbits
% planet_id  - integer from 1 to 7 (from Mercury to Uranus)
% date       - [year,month,day,hour,minute,second]
% col        - color spec
% linewidth  - linewidth spec
%
% -----------------------------------------------------------------------%


muSun=132712439935; % sun's gravitational parameter
[~,R,V,~]=planet_elements_and_sv(assumption,planet_id,date(1),date(2),date(3),date(4),date(5),date(6) );
OrbPlot(R,V,muSun,col, linewidth)
