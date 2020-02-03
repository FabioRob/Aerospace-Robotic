function [coe , rd,vd , ra,va ,v_planet_dep , v_planet_arr ,v_inf_dep ,v_inf_arr  ] = lambert_trajectory ( planet_id_dep , planet_id_arr , date_dep , date_arr ,plotOrbit ,linewidth)

% -----------------------------------------------------------------------%
%
% lambert_trajectory function resolves lambert problem from two planet given the date of departure
% and date of arrival, with assumption of coplanar orbits.
% Frame is heliocentric ecliptic frame.
%
% Arguments :
%
% planet_id_dep  - departure planet,from 1 to 7 (from Mercury to Uranus)
% planet_id_arr  - arrival planet, from 1 to 7 (from Mercury to Uranus)
% date_dep       - departure date([year,month,day,hour,minute,second])
% date_dep       - arrival date([year,month,day,hour,minute,second])
% plotOrbit      - if it's true, it plots the trajectory
% linewidth      - linewidth spec
%
% Output :
%
% coe            - orbital elements at arrival planet [h e RA incl w TA a ] with all angles in radians
% rd,vd          - state vector at departure 
% ra,va          - state vector on arrival
% v_planet_dep   - velocity vector of departure planet
% v_planet_arr   - velocity vector of arrival planet
% v_inf_dep      - relative velocity at departure
% v_inf_arr      - relative velocity on arrival
%
% -----------------------------------------------------------------------%

muSun=132712439935; % sun's gravitational parameter
deg = pi/180;

[~,rd,v_planet_dep, jd_dep ] = planet_elements_and_sv( 3, planet_id_dep , date_dep(1),date_dep(2),date_dep(3),date_dep(4),date_dep(5),date_dep(6)) ;
[coe_pl_arr ,ra,v_planet_arr, jd_arr] = planet_elements_and_sv(3,planet_id_arr, date_arr(1),date_arr(2),date_arr(3),date_arr(4),date_arr(5),date_arr(6));

tof=(jd_arr-jd_dep)*86400;

[vd,va] = lambert(rd,ra, tof, 'pro') ; 

coe = coe_from_sv(ra,va,muSun);
TA_arr=(zero_to_360( coe(6)/deg ) )*deg ; % true anomaly at arrival(rad)
w= (coe_pl_arr(5) + coe_pl_arr(6) - TA_arr/deg) *deg; % argument of perihelium(rad)
coe(5) = w ;
a = coe(7); % semimajor axis

v_inf_dep = vd-v_planet_dep ; % Excess velocity vector at departure
v_inf_arr =  va-v_planet_arr ; % Excess velocity vector at arrival

if plotOrbit==true
    OrbPlot(rd,vd,muSun,'r',linewidth,TA_arr)
    %plotApseLine(coe,muSun,'r')
end



