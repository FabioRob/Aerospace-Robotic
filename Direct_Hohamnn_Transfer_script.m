% This script computes total deltaV and time of flight, and plots all figures
% (3 figures : departure fase,cruise fase and arrival fase)
% of a Direct Hohmann Transfer from Earth to Uranus, with an initial 200 Km
% altitude parking orbit around Earth and an arrival 1000 Km orbit around
% Uranus.

initial_date=[2020,1,1,0,0,0]; % launch must occurs after this date(progect specification 2)

[muEarth,radius_E,~] = physical_data(3) ; 
[muUranus,radius_U,~] = physical_data(7) ;
deg=pi/180; % to convert from degrees to radians

% Compute the Hohmann transfer
[ v_inf_dep , v_dep , v_inf_arr , v_arr , h_H , e_H , a_H , TOF_H  ] = Compute_Direct_Hohmann ( 3 , 7  ) ;

% iteratevely search for a date(starting from 1st January 2020)
% in which phase angle between the two
% planets is such that the spacecraft correctly rendezvous on Uranus
% phase must be equal to phase_in
[date_dep,date_arr ] = Find_Direct_Hohmann_Departure_Date ( 3 , 7 , initial_date ) ;

[V_E_H , V_U_H ] = plot_Hohmann_Transfer(date_dep , date_arr) ;


% Departure from circular ,200 km altitude, equatorial, Earth's parking  orbit
z_E_P = 200; % parking orbit altitude
r_E_P = radius_E + z_E_P; % parking orbit radius
v_E_P = sqrt(muEarth / r_E_P); % velocity on parking orbit

% First of all, rotate parking orbit plane to the ecliptic plane

E_obl = 23.45 * deg ; % Earth's equatorial plane angle to ecliptic (rad)
dV_plane_change_E = 2 * v_E_P * sin(E_obl/2) ; % deltaV needed for the maneuver

[ dV_hyp_dep , TOF_exit ] = EarthDeparture(z_E_P, V_E_H , date_dep , 2 , true );

%%%%%%%%%%%%%%%%%%%%%%%%%

% Arrival on Uranus' orbit at 1000 km altitude
z_U_P = 1000 ; %  orbit altitude
r_U_P = radius_U + z_U_P ; %  orbit radius
v_U_P = sqrt(muUranus / r_U_P) ; % velocity in the circular orbit

[dV_hyp_arr , TOF_in ] = UranusCapture ( z_U_P , 0 , V_U_H , date_arr , 2 , true ) ;

% Finally, rotate orbit plane to polar orbit around Uranus

U_obl = 97.77 * deg ; % Uranus's equatorial plane angle to ecliptic (rad)
U_incl =abs( pi/2 - U_obl ) ; % Inclination between ecliptic plane an the final polar orbit (rad)
dV_plane_change_U = 2 * v_U_P * sin(U_incl/2) ; % deltaV needed for the maneuver


dV_plane_changes = dV_plane_change_E + dV_plane_change_U ; % total dV for plane changes at Eartj and Uranus

dV_HohmannTransfer = dV_hyp_dep + dV_hyp_arr ; % total dV for the direct Hohmann transfer from a elciptic 200 km circular Earth orbit 
 % to a ecliptic 1000 km circular Uranus orbit

deltaV_Hohmann = dV_plane_changes + dV_HohmannTransfer ; % total dV using a direct Hohmann transfer(km/s)

TOF_Hohmann = (TOF_exit + TOF_H + TOF_in) / 86400 / 365 ; % total time of flight(years)

disp('Total deltaV needed for the Hohmann Transfer is(Km/s): ')
deltaV_Hohmann
disp('Total approximated time of the Hohmann Transfer is(years) : ')
TOF_Hohmann
