% This script computes total deltaV and time of flight of the mission to Uranus with
% three flybys VEE(Venus-Earth-Earth) , and plots all figures :
% figure 1 : heliocentric planets and spacecraft orbits 
% figure 2 : trajectory from second Earth FlyBy to Uranus
% figure 3 : departure form Earth
% figure 4 : Venus Flyby
% figure 5 : First Earth Flyby
% figure 6 : Second Earth Flyby
% figure 7 : Uranus elliptical capture orbit
% figure 8 : Uranus circular capture orbit(greater dV needed)

load 'planets_data'

figure
hold all
title('VEE FlyBys to Uranus')
whitebg('k')
axis equal
scatter3(0,0,0,200,'y','filled') % Sun

date_dep = [2020 , 3 , 1 ,0,0,0] ; % departure from Earth

jd_dep = date2JD( date_dep(1),date_dep(2),date_dep(3),date_dep(4),date_dep(5),date_dep(6));
plotPlanetOrbit (3 , 3 , date_dep , 'w' , 1 )
plotPlanetPosition (3 ,3 ,date_dep , 'o' , 8 , 'w' )
z_parking = 200; % equatorial parking orbit altitude
v_E_P=sqrt(muEarth/(radius_E+z_parking));

plotPlanetOrbit (3 , 4 , date_dep , 'r' , 1 ) % Mars orbit

%  Earth-->Venus and Venus flyby

date_fb_Venus = [2020 , 9 ,20 , 0 ,0 ,0] ; % Venus flyby date


jd_fb_Venus = date2JD( date_fb_Venus(1),date_fb_Venus(2),date_fb_Venus(3),date_fb_Venus(4),date_fb_Venus(5),date_fb_Venus(6));
tof_E_V = (jd_fb_Venus - jd_dep ) *86400 ;
plotPlanetOrbit (3 , 2 , date_fb_Venus , 'g' , 2 )
view([0 90])
plotPlanetPosition(3, 2 ,date_fb_Venus , '*' , 8 , 'w' )
[coe_E_V , rd_s_E ,vd_s_E , ra_s_V ,va_s_V ,v_E_dep , v_V_arr ,v_inf_dep_E ,v_inf_arr_V  ] = lambert_trajectory ( 3 , 2 , date_dep , date_fb_Venus ,  true , 2) ;

coe_dep = coe_from_sv(rd_s_E,vd_s_E,muSun);
TA_dep = coe_dep(6); 

E_obl = 23.45 * deg ; % Earth's equatorial plane angle to ecliptic (rad)
dV_plane_change_E = 2 * v_E_P * sin(E_obl/2) ; % deltaV needed for the maneuver

% Find altitude of flyby and date which will make the spacecraft encounter Earth to do a Earth's flyby
[date_fb_Earth_1 , zp_fb_V ] = fb_date_and_altitude ( 2 , 3 , date_fb_Venus , [2021 , 3 ,1 , 0 ,0 ,0], [2021 , 12 ,1 , 0 ,0 ,0] , 0 , 1e-4 , coe_E_V ) ;

[ dVgained_1 , coe_out_1 , v_out_1 ,TOF_fb_1] = FlyBy ( 2 , date_fb_Venus , coe_E_V  , zp_fb_V ,false,1) ;
[~, ra_V, va_V, ~ ] = planet_elements_and_sv(3,2, date_fb_Venus(1),date_fb_Venus(2),date_fb_Venus(3),date_fb_Venus(4),date_fb_Venus(5),date_fb_Venus(6));

jd_fb_Earth_1 = date2JD( date_fb_Earth_1(1),date_fb_Earth_1(2),date_fb_Earth_1(3),date_fb_Earth_1(4),date_fb_Earth_1(5),date_fb_Earth_1(6));
tof_V_E = (jd_fb_Earth_1 - jd_fb_Venus ) * 86400 ; % sec
[r_s_arr_Earth_1 , v_s_arr_Earth_1 ] = rv_from_r0v0(ra_s_V, v_out_1 , tof_V_E , muSun) ;
plot3(r_s_arr_Earth_1(1),r_s_arr_Earth_1(2),r_s_arr_Earth_1(3),'*','markerEdgeColor','w', 'markersize',10)
plotPlanetPosition(3, 3 ,date_fb_Earth_1 , 'o' ,10, 'w' )
[~, ra_E_1, va_E_1, ~ ] = planet_elements_and_sv(3,3, date_fb_Earth_1(1),date_fb_Earth_1(2),date_fb_Earth_1(3),date_fb_Earth_1(4),date_fb_Earth_1(5),date_fb_Earth_1(6));
distance_spacecraft_Earth_1 = norm(r_s_arr_Earth_1 - ra_E_1 ) ;  % it must be less then SOI_E = 925000 km ( 0.925e6 )

coe_in_E_1 = coe_from_sv(r_s_arr_Earth_1,v_s_arr_Earth_1,muSun);

OrbPlot(ra_s_V,v_out_1 ,muSun ,'b' ,3 , coe_in_E_1(6)) % post Venus flyby orbit

% SECOND FLYBY
% first Earth flyby to get to Earth again

% Find altitude of flyby and date which will make the spacecraft encounter Earth again to do a second Earth's flyby
[date_fb_Earth_2 , zp_fb_E_1 ] = fb_date_and_altitude ( 3 , 3 , date_fb_Earth_1 , [2024 , 1 ,1 , 0 ,0 ,0] , [2026 , 1 ,1 , 0 ,0 ,0] , 0 , 1e-5  , coe_in_E_1 ) ;

[ dVgained_2 , coe_out_2 , v_out_2 ,TOF_fb_2 ] = FlyBy ( 3 , date_fb_Earth_1 , coe_in_E_1  , zp_fb_E_1 , false,2 ) ; % First Earth flyby

T_out_2 =  ((2*pi/sqrt(muSun))*(coe_out_2(7))^1.5 ) / 86400 / 365  ; % post flyby orbit's period in years

OrbPlot(ra_E_1,v_out_2 ,muSun ,'y' ,3 ) % post first Earth-flyby orbit

jd_fb_Earth_2 = date2JD( date_fb_Earth_2(1),date_fb_Earth_2(2),date_fb_Earth_2(3),date_fb_Earth_2(4),date_fb_Earth_2(5),date_fb_Earth_2(6));
plotPlanetPosition(3, 3 ,date_fb_Earth_2 , 'o' , 3 , 'w' )

tof_E_E = (jd_fb_Earth_2 - jd_fb_Earth_1 ) * 86400 ; % sec
[r_s_arr_Earth_2 , v_s_arr_Earth_2 ] = rv_from_r0v0(ra_E_1, v_out_2 , tof_E_E , muSun) ;
plot3(r_s_arr_Earth_2(1),r_s_arr_Earth_2(2),r_s_arr_Earth_2(3),'d','markerEdgeColor','w', 'markersize',3)
[~, ra_E_2, va_E_2, ~ ] = planet_elements_and_sv(3,3, date_fb_Earth_2(1),date_fb_Earth_2(2),date_fb_Earth_2(3),date_fb_Earth_2(4),date_fb_Earth_2(5),date_fb_Earth_2(6));
distance_spacecraft_Earth_2 = norm(r_s_arr_Earth_2 - ra_E_2 ) ;  % it must be less then SOI_E = 925000 km ( 0.925e6 )

coe_in_E_2 = coe_from_sv(r_s_arr_Earth_2,v_s_arr_Earth_2,muSun);

% Third Flyby
% Second Earth flyby

date_Uranus=[2038,6,1,0,0,0] ;
zp_fb_E_2 = 200 ; 

[ dV_gained_3 , coe_out_3 , v_out_3 ,TOF_fb_3 ] = FlyBy ( 3 , date_fb_Earth_2 , coe_in_E_2  , zp_fb_E_2 , false , 3 ) ; % Second Earth FlyBy
T_out_3 =  ((2*pi/sqrt(muSun))*(coe_out_3(7))^1.5 ) / 86400 / 365 ;
OrbPlot(ra_E_2,v_out_3 ,muSun ,'b' , 1 ,1.5  ) % post second Earth-flyby orbit(not travelled)

lambert_trajectory ( 3 , 7 , date_fb_Earth_2 , date_Uranus ,  true , 2) ;
text(2e8,2e8,'to Uranus\rightarrow')


xlim([-5e8,3e8])
ylim([-2e8,4e8])

text(-9e7,-3e7,'\leftarrowVenus orbit')

text(-15e7,14e7,'Earth orbit\rightarrow')

text(-28e7,5e7,[ 'Launch (',num2str(date_dep(3)),'/',num2str(date_dep(2)),'/',num2str(date_dep(1)),')\rightarrow'])

text(4e7,9e7,'\uparrow')
text(-1e7,8e7,'Venus FlyBy')
text(-1e7,7e7,[ '(',num2str(date_fb_Venus(3)),'/',num2str(date_fb_Venus(2)),'/',num2str(date_fb_Venus(1)),')'])

text(-9e7,-9e7,'E2V')
text(-8e7,-10e7,'\downarrow')

text(-3.3e8,-5e7,'Post Venus FlyBy orbit\rightarrow')

text(9e7,-13e7,['\leftarrow1_{st}Earth FlyBy(',num2str(date_fb_Earth_1(3)),'/',num2str(date_fb_Earth_1(2)),'/',num2str(date_fb_Earth_1(1)),')'])
text(-4e8,2.5e8,'\leftarrowPost 1_{st}Earth FlyBy orbit')

text(7e7,-15e7,'\uparrow')
text(5e7,-16e7,['2_{nd}Earth FlyBy(',num2str(date_fb_Earth_2(3)),'/',num2str(date_fb_Earth_2(2)),'/',num2str(date_fb_Earth_2(1)),')'])

%%
% to Uranus
figure
hold all
axis equal
whitebg('k')


vel_scale = 5e7; % scale factor for velocity vectors

title('Heliocentric trajectory to Uranus after second Earth FlyBy ')

scatter3(0,0,0,200,'y','filled') % Sun
plotPlanetOrbit (3 , 3 , date_fb_Earth_2 , 'w' , 1 )
plotPlanetPosition (3 ,3 ,date_fb_Earth_2 , 'o' , 6 , 'w' )

[~ , rd_s ,vd_s , ~ ,~ ,~ , ~ ,~ ,~  ] = lambert_trajectory ( 3 , 7 , date_fb_Earth_2 , date_Uranus ,  true , 2) ;

xlim([-2e9,2e9])
ylim([-0.5e9,3e9])

plotPlanetOrbit (3 , 4 , date_dep , 'r' , 1 ) % Mars orbit
plotPlanetOrbit (3 , 5 , date_dep , 'b' , 1 ) % Jupiter orbit
plotPlanetOrbit (3 , 6 , date_dep , 'y' , 1 ) % Saturn orbit

[~, ra_U, v_U_arr, jd_Uranus ] = planet_elements_and_sv(3,7, date_Uranus(1),date_Uranus(2),date_Uranus(3),date_Uranus(4),date_Uranus(5),date_Uranus(6));

deltaV_burn  = deltaV_target (3 ,date_fb_Earth_2 , coe_out_3 , v_out_3 , 7 , date_Uranus ) ;

tof_E_U = (jd_Uranus - jd_fb_Earth_2)*86400 ;

[~,va_s_U] = rv_from_r0v0( rd_s ,vd_s , tof_E_U , muSun) ;

plotPlanetOrbit (3 , 7 , date_Uranus , 'c' , 2 )
plotPlanetPosition (3 , 7 ,date_Uranus , 'o' , 10 , 'w' )

quiver(ra_E_2(1),ra_E_2(2),v_out_3(1)*vel_scale,v_out_3(2)*vel_scale,'color','b')
text(5e8,7e8,'Velocity after 2_{nd}Earth FlyBy')
quiver(ra_E_2(1),ra_E_2(2),vd_s(1)*vel_scale,vd_s(2)*vel_scale,'color','r')
text(8e8,2e8,'Velocity in final orbit to Uranus')
quiver(ra_E_2(1)+0.9*v_out_3(1)*vel_scale,ra_E_2(2)+0.9*v_out_3(2)*vel_scale,(vd_s(1)-v_out_3(1))*vel_scale,(vd_s(2)-v_out_3(2))*vel_scale,'color','g')
text(15e8,8e8,[ 'dV(',num2str(deltaV_burn,2), 'Km/s )' ] )

text(-5e8,0,'Earth orbit\rightarrow')
text(1e9,26e8,'\uparrow')
text(8e8,25e8,'Uranus orbit')
text(5e7,-20e7,['2_{nd}Earth FlyBy(',num2str(date_fb_Earth_2(3)),'/',num2str(date_fb_Earth_2(2)),'/',num2str(date_fb_Earth_2(1)),')'])

text(-12e8,27e8,['Arrival on Uranus(',num2str(date_Uranus(3)),'/',num2str(date_Uranus(2)),'/',num2str(date_Uranus(1)),')'])

% Uranus Capture(elliptical orbit)
%reduce deltaV of capture using an ellipitical captur orbit

v_inf_arr = norm(va_s_U - v_U_arr ) ;

% Arrival on Uranus' orbit at 1000 km altitude
z_U_P = 1000 ; %  final orbit altitude

% Instead of capturing directly into a circular orbit, capture orbit has a
% eccentricity of 0.5, and then an impulsive maneuver is performed to
% change eccentricity to zero(see p. 225 Mengali)
% In this way, total deltaV of arrival is obviously the same, but it's splitted
% into two differents impulsive maneuvers so it could be more realistic to
% achieve.

e_capture = 0.5; % eccentricity of capture orbit
r_U_P = radius_U + z_U_P ; % capture orbit periapsis radius

[deltaV_arrival_elliptic , ~ ] = UranusCapture ( z_U_P , e_capture , va_s_U , date_Uranus , 3 , false ) ;

r_U_A = r_U_P * (1+e_capture)/(1-e_capture); % apoapsis radii of elliptic capture orbit
v_U_P_capt = sqrt(2*( muUranus/r_U_P - muUranus/(r_U_A+r_U_P) ) ); % capture orbit's periapsis velocity
v_U_P_circ = sqrt(muUranus / r_U_P) ; % velocity in a circular orbit with periapsis altitude of 1000 km

dV_to_circ = v_U_P_capt - v_U_P_circ ; % second impulsive manuever to change from elliptical to circular orbit

dV_arr_tot = deltaV_arrival_elliptic + dV_to_circ ; % equal to deltaV for a circular captur orbit 1000Km altitude


%  Uranus Capture(circular orbit)

z_U_P = 1000 ;
r_U_P = radius_U + z_U_P ; % capture orbit periapsis radius
v_U_P = sqrt(muUranus / r_U_P);

% compute deltaV of departure,arrival and plot the orbits
[deltaV_departure , TOF_exit ] = EarthDeparture(z_parking, vd_s_E , date_dep , 3 ,false ) ;
[deltaV_arrival , TOF_in ] = UranusCapture ( z_U_P , 0 , va_s_U , date_Uranus , 3 , false ) ;

% Finally, rotate orbit plane to polar orbit around Uranus

U_obl = 97.77 * deg ; % Uranus's equatorial plane angle to ecliptic (rad)
U_incl =abs( pi/2 - U_obl ) ; % Inclination between ecliptic plane an the final polar orbit (rad)
dV_plane_change_U = 2 * v_U_P * sin(U_incl/2) ; % deltaV needed for the maneuver

deltaV = dV_plane_change_E + deltaV_departure + deltaV_burn + deltaV_arrival + dV_plane_change_U ;

TOF = ( tof_E_V + tof_V_E + tof_E_E + tof_E_U ) / 86400 / 365 ;  % years

%TOF_tot = (TOF + TOF_exit + TOF_fb_1 +TOF_fb_2 + TOF_fb_3 + TOF_in)  / 86400 / 365  % years

%%
% Create planetary departure,arrival, and flybys figures
 EarthDeparture(z_parking, vd_s_E , date_dep , 3 , true ) ; % Earth departure
 
 FlyBy ( 2 , date_fb_Venus , coe_E_V  , zp_fb_V ,true,1) ; % Venus FlyBy
 FlyBy ( 3 , date_fb_Earth_1 , coe_in_E_1  , zp_fb_E_1 , true,2 ) ; % First Earth FlyBy
 FlyBy ( 3 , date_fb_Earth_2 , coe_in_E_2  , zp_fb_E_2 , true , 3 ) ; % Second Earth FlyBy

 UranusCapture ( z_U_P , e_capture , va_s_U , date_Uranus , 3 , true ) ; % Uranus capture in an elliptic orbit
 UranusCapture ( z_U_P , 0 , va_s_U , date_Uranus , 3 , true ) ; % Uranus capture in a circular orbit

 disp('Total deltaV needed for this mission is(Km/s): ')
 deltaV
 disp('Total approximated time of flight of this mission is(years) : ')
 TOF
