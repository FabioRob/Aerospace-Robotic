function [deltaV_arr , TOF_in ] = UranusCapture ( z_capture , e_capture , V_arr , date , assumption , graphic )

% -----------------------------------------------------------------------%
%
% UranusCapture function computes deltaV burn needed and TOF of Uranus capture 
%
% Arguments:
%
% z_capture  - periapsis altitude of capture orbit
% e_capture  - eccentricity of capture orbit
% V_arr      - heliocentric velocity of the spacecraft at arrival on Uranus
% date       - arrival date([year,month,day,hour,minute,second])
% assumption - 1/2/3 for,respectively, assumption of real/coplanar&circular/coplanar orbits
% graphic    - boolean,if it's true, it plots arrival hyperbola
%
% Output:
%
% deltaV_arr - deltaV burn necessary to enter in Uranus capture orbit
% TOF_in     - time of flight of arrival hyperbola
%
% -----------------------------------------------------------------------%


% Earth's gravitational parameter(Km^3/s^2) , Earth's radius(Km) and
% Earth's SOI radius(Km)
[muUranus,radius_U,SOI_U] = physical_data(7) ;
deg=pi/180; % to convert from degrees to radians

[~, ~, v_U_arr, ~] = planet_elements_and_sv(assumption,7, date(1) ,date(2), date(3), date(4), date(5), date(6));
v_U = norm(v_U_arr); % heliocentric speed of Uranus at arrival date

v_inf_arr = norm(V_arr - v_U_arr ) ;

r_U_P = radius_U + z_capture ; %  capture orbit periapsis radius
v_U_P = sqrt(muUranus*(1+e_capture) / r_U_P) ; % velocity at periapsis of capture orbit

% Arrival hyperbola

h_hyp_arr = r_U_P * sqrt ( v_inf_arr^2 + 2 *muUranus/r_U_P )  ; % specific angular momentum
e_hyp_arr = 1 + r_U_P*v_inf_arr^2/muUranus ; % eccentricity
TA_inf_arr = acos(-1/e_hyp_arr) ; % asyptotic value of the true anomaly
beta_hyp_arr = pi - TA_inf_arr ; % angle between asymptote and apse line
a_hyp_arr = h_hyp_arr^2 / (muUranus * (e_hyp_arr^2-1)) ;  % semimajor axis
aiming_hyp_arr = a_hyp_arr * sqrt(e_hyp_arr^2 -1 )  ; % aiming radius
v_hyp_per_arr = h_hyp_arr / r_U_P ;  % velocity at perigeo

% Delta V needed for Uranus' SOI capture
deltaV_arr = v_hyp_per_arr - v_U_P ;

% Time of flight from Uranus SOI entrance to hyperbola's periapsis 

TA_in = -acos ( e_hyp_arr^-1 * (h_hyp_arr^2/(muUranus*SOI_U) -1 ) ) ; % true anomaly at SOI entrance
TA_exit = 0 ;  % true anomaly at arrival on parking orbit
F_in = etheta2E(e_hyp_arr, TA_in); % eccentric anomaly at SOI entrance
M_in = e_hyp_arr * sinh(F_in) - F_in ; % mean anomaly at SOI entrance
n_in = sqrt ( muUranus / a_hyp_arr^3 ) ; % mean motion at SOI entrance
TOF_in = abs(M_in / n_in) ;

if graphic==true
    % Plot Uranus, arrival orbit on the ecliptic and the capture hyperbola
    % XY frame is a geocentric frame with Y directed as the planet orbit track
    % X directed in opposite direction of the sun(vernal equinox line)
    
    figure
    if e_capture==0
        title(['Uranus Capture(circular capture orbit ',num2str(date(3)),'/',num2str(date(2)),'/',num2str(date(1)),')'])
    else
        title(['Uranus Capture(elliptical capture orbit ',num2str(date(3)),'/',num2str(date(2)),'/',num2str(date(1)),')'])
    end
    
    hold all
    whitebg('k') 
    scatter3(0,0,0,30,'w','filled')

    % Compute rotation matrix from perifocal frame of the hyperbola(and capture orbit)
    % to planetocentric frame
    RA_hyp_arr = 0; 
    i_hyp_arr  = 0; 
    w_hyp_arr = (pi/2 - beta_hyp_arr) ; % argument of periapsis(from X axis to perifocal x axis)
    R_from_per_arr  =  (R_to_perifocal(RA_hyp_arr ,i_hyp_arr ,w_hyp_arr ))' ; % rotation matrix from peifocal frame to planetocentric frame

    R_p_arr  = R_from_per_arr  * [ r_U_P , 0 ,0 ]' ; % position vector of the spacecraft at periapsis, in geocentric frame
    V_p_arr  = R_from_per_arr *[0,v_hyp_per_arr,0]' ;  % Velocity vector of the spacecraft at periapsis' hyperbola in g. frame
    V_p_capt  = R_from_per_arr *[0,v_U_P,0]' ;  % Velocity vector of the spacecraft in capture orbit
    OrbPlot( R_p_arr  , V_p_arr , muUranus ,'g' , 1 , TA_in ) % plot hyperbola
    
    OrbPlot ( R_p_arr , V_p_capt , muUranus , 'r' , 1 ) % plot capture orbit
    text(-1e5,0,'Capture orbit 1000km')
    text(-8e4,-0.5e4,'altitude')

    
    xlim([-1e5  1e5 ])
    ylim([-1e5  1e5 ])
    axis equal
    
    plot3([0;R_p_arr(1)],[0;R_p_arr(2)],[0;R_p_arr(3)],'--') % draw periapsis radius of hyperbola
    scatter3(R_p_arr(1),R_p_arr(2),R_p_arr(3),50,'w','filled') % draw point in which the maneuver occurs

    vel_scale_U = 6e3 ; % scale factor to visualize velocity vectors at Uranus capture
    %viscircles([0,0] ,SOI_U) 

    quiver3(3e4,0,0, 3e4,0,0, 'Color','y','LineWidth',1) % direction to the Sun
    text(3e4,1e4,'to the Sun')

    
    quiver3(0,0,0,0,-v_U*vel_scale_U,0,'Color','b','LineWidth',1) % Uranus velocity
    text(-5e4,-5e4,['Uranus velocity=',num2str(v_U), ' km/s'])

    quiver3(aiming_hyp_arr,-1e5,0 , 0,v_inf_arr*vel_scale_U,0 ,'Color' ,'g' , 'LineWidth' ,1) % Excess speed
    text(6e4,-6e4,['Excess speed=',num2str(v_inf_arr) ,' km/s'])

    ang([0,0],8e3,[0 pi/2-beta_hyp_arr],'r') % angle between hyperbola's apse line and normal to the asympote
    text(2e3,0,['\pi/2-\beta=',num2str(w_hyp_arr/deg),'\circ'])

    plot3([aiming_hyp_arr;aiming_hyp_arr],[-1e5;1.5e5],[0;0],'g--') % asymptote
    text(aiming_hyp_arr-5e4,6e4,'asymptote')

    quiver3(R_p_arr(1),R_p_arr(2),R_p_arr(3),V_p_capt(1)*vel_scale_U,V_p_capt(2)*vel_scale_U,V_p_capt(3)*vel_scale_U,'Color','r','LineWidth' ,1)
    text(-2e4,6e4,'spacecraft velocity')
    text(-2e4,5e4 ,'in capture orbit')
    text(-2e4,4e4,['(',num2str(v_U_P) ,' km/s )'])


    quiver3(R_p_arr(1)+1e3,R_p_arr(2)-1e3*tan(pi/2-beta_hyp_arr),R_p_arr(3),V_p_arr(1)*vel_scale_U,V_p_arr(2)*vel_scale_U,V_p_arr(3)*vel_scale_U,'Color','g','LineWidth' ,1)
    text(-1e5,4e4,'spacecraft velocity')
    text(-1e5,3e4,'in hyperbolic orbit')
    text(-1e5,2e4,['(',num2str(v_hyp_per_arr) ,' km/s )'])

    text(-5e4,8e4,['\fontsize{15}\DeltaV capture = ',num2str(deltaV_arr),' km/s'])
    
end


