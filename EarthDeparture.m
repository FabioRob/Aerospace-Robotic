function [ deltaV , TOF_exit ] = EarthDeparture(z_parking, V_dep , date , assumption ,graphic )

% -----------------------------------------------------------------------%
%
% EarthDeparture function computes hyperbola of departure from circular orbit around Earth
%
% Arguments:
%
% z_parking  - altitude of parking orbit(assumed circular)
% V_dep      - heliocentric velocity vector at departure from Earth(row vector)
% date       - departure date from Earth([year,month,day,hour,minute,second])
% assumption - 1,2,3 for,repectively, real,coplanar and circular, coplanar orbits
% graphic    - boolean, if it's true, it plots parking orbit, hyperbola of departure and velocity vectors
%
% Outputs:
%
% deltaV     - velocity burn needed to put the spacecraft in the hyperbolic
%              departure trajectory
% TOF_exit   - time of flight from huperbola's periapsis to the exit from
%              Earth's SOI
%
% -----------------------------------------------------------------------%

% Earth's gravitational parameter(Km^3/s^2) , Earth's radius(Km) and
% Earth's SOI radius(Km)
[muEarth,radius_E,SOI_E] = physical_data(3) ; 
deg = pi/180 ; % to convert from degrees to radians

% heliocentric velocity vector of Earth
[~, ~, V_Earth, ~] = planet_elements_and_sv(assumption,3, date(1) ,date(2), date(3), date(4), date(5), date(6));
v_E = norm(V_Earth); % heliocentric speed of Earth at departure date
r_p = radius_E + z_parking ; % parking orbit radius(equal to hyperbola's periapsis radius) 
v_inf = norm(V_dep-V_Earth ) ; % excess speed 
v_parking = sqrt(muEarth/r_p) ; % velocity in parking orbit

e = 1 + r_p * v_inf^2 / muEarth ; % departure hyperbola's eccentricity
h = r_p * sqrt( v_inf^2 + 2*muEarth / r_p ) ; % departure hyperbola's specific angular momentum (km^2/s)

v_p = h / r_p ; % velocity at hyperbola's periapsis

deltaV = v_p - v_parking ; % speed burn required

TA_inf = acos(-1/e) ; % asymptotic value of the true anomaly
beta = pi - TA_inf ; % angle between asymptote and apse line
delta = 2*asin(1/e); % turn angle
a = h^2 / (muEarth * (e^2-1)) ;  % semimajor axis
aiming = a * sqrt(e^2 -1 )  ; % aiming radius

% Time of flight from hyperbola's periapsis to exit from Earth's SOI
TA_exit = acos ( e^-1 * (h^2/(muEarth*SOI_E) -1 ) ) ;  % true anomaly at SOI exit
F_exit = etheta2E(e, TA_exit); % eccentric anomaly at SOI exit
M_exit = e * sinh(F_exit) - F_exit ; % mean anomaly at SOI exit
n_exit = sqrt ( muEarth / a^3 ) ; % mean motion at SOI exit
TOF_exit = M_exit / n_exit ; % Time of flight to exit SOI(sec)

% hyperbola's Right Ascension,inclination and argument of perigee
% (assuming that the hyperbolic escape trajectory lies in the same plane of
% the parking orbit)
RA = 0; 
i = 0; 
w= -(pi/2 - beta) ; % argument of perigeo(from X axis to hyperbola's perifocal x axis)
R_from_per =  (R_to_perifocal(RA,i,w))' ; % rotation matrix from hyperbola' perifocal frame to geocentric frame

R_p = R_from_per * [ r_p , 0 ,0 ]' ; % position vector of the spacecraft at periapsis, in geocentric frame
V_p = R_from_per*[0,v_p,0]' ;  % Velocity vector of the spacecraft at periapsis' hyperbola in g. frame
V_parking = R_from_per*[0,v_parking,0]' ;  % Velocity vector of the spacecraft in circular orbit

if graphic==true
    % Create Earth Departure plot
    figure
    title(['Earth Departure (',num2str(date(3)),'/',num2str(date(2)),'/',num2str(date(1)),')'])
    hold all
    whitebg('k')    
    scatter3(0,0,0,30,'w','filled')

    OrbPlot ( [r_p,0,0] , [0,v_parking,0] , muEarth , 'r' , 1) % plot spacecraft's parking orbit
    text(-8e3,-8e3,['Circular orbit ',num2str(z_parking), ' km altitude'])

    OrbPlot( R_p , V_p, muEarth ,'g' , 1 ,  TA_exit  ) % plot escape hyperbola

    xlim([-2e4  2e4 ])
    ylim([-2e4  2e4 ])
    axis equal

    plot3([0;R_p(1)],[0;R_p(2)],[0;R_p(3)],'--') % draw periapsis radius of hyperbola
    scatter3(R_p(1),R_p(2),R_p(3),50,'w','filled') % draw point in which the maneuver occurs

    vel_scale_E = 0.5e3 ; % scale factor to visualize velocity vectors at Earth departure     

    quiver3(-0.3e4,0,0, -0.5e4,0,0, 'Color','y','LineWidth',1) % direction to the Sun
    text(-0.6e4,5e2,'to the Sun')

    quiver3(0,0,0,0,v_E*vel_scale_E,0,'Color','b','LineWidth',1) % Earth velocity
    text(-5e3,14e3,['Earth velocity=',num2str(v_E), ' km/s'])

    quiver3(aiming,8e3,0 , 0,v_inf*vel_scale_E,0 ,'Color' ,'g' , 'LineWidth' ,1) % Excess speed
    text(1e4,1e4,['Excess speed=',num2str(v_inf) ,' km/s'])

    an = (pi/2 - beta)/deg;
    ang([0,0],2e3,[0 -pi/2+beta],'r') % angle between hyperbola's apse line and normal to the asympote
    text(2e3,0,['\pi/2-\beta=',num2str(an),'\circ'])

    plot3([aiming;aiming],[-1e4;1.5e4],[0;0],'g--') % asymptote
    text(aiming,-8e3,'asymptote')

    % Spacecraft velocity in circular orbit
    quiver3(R_p(1),R_p(2),R_p(3),V_parking(1)*vel_scale_E,V_parking(2)*vel_scale_E,V_parking(3)*vel_scale_E,'Color','r','LineWidth' ,1)
    text(2e3,-3e3,'spacecraft velocity')
    text(2e3,-4e3 ,'in circular orbit')
    text(2e3,-5e3,['(',num2str(v_parking) ,' km/s )'])

    % Spacecraft velocity in hyperbolc orbit at periapsis
    quiver3(R_p(1)+1e3,R_p(2)-1e3,R_p(3),V_p(1)*vel_scale_E,V_p(2)*vel_scale_E,V_p(3)*vel_scale_E,'Color','g','LineWidth' ,1)
    text(8e3,1e3,'spacecraft velocity')
    text(8e3,0,'in hyperbolic orbit')
    text(8e3,-1e3,['(',num2str(v_p) ,' km/s )'])

    % DeltaV needed for the maneuver
    text(-1.5e4,1e4,['\fontsize{15}\DeltaV departure = ',num2str(deltaV),' km/s'])

end