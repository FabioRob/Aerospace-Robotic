function [ dVgained , coe_out , v_out_xy , TOF_fb ] = FlyBy ( planet_id , date , coe_in  , zp ,graphic,fb_number)

% This function computes flyby trajectory
%
% Arguments:
%
% planet_id - from 1 to 7 (from Mercury to Uranus)  
% date      - flyby date([year,month,day,hour,minute,second])
% coe_in    - orbital elements at arrival ( [h e RA incl w TA a ] with all angles in radians )
% zp        - flyby hyperbola's periapsis altitude from planet
% graphic   - boolean, if it's true, it plots flyby hyperbola
% fb_number - 1/2/3 to plot,respectively, Venus/Earth/Second Earth flyby hyperbolas
%
% Outputs:
%
% dVgained  - deltaV gained through the flyby
% coe_out   - orbital elements of the post flyby orbit
% v_out_xy  - heliocentric velocity(in heliocentric frame) at exit from SOI
% TOF_fb    - flyby's time of flight
%
% -----------------------------------------------------------------------%


deg = pi/180 ; % to convert from degrees to radians
muSun=132712439935; % sun's gravitational parameter
[mu,radius,SOI] = physical_data(planet_id);

% periapsis of the flyby hyperbola
rp = radius + zp ;

% Orbital elements and state vector of the planet
[coe_pl , R, V , ~] = planet_elements_and_sv(3,planet_id, date(1),date(2),date(3),date(4),date(5),date(6) ) ;
% Path angle of planet velocity
path_angle_planet = atan2(coe_pl(2)*sin(coe_pl(6)*pi/180) , 1 + coe_pl(2)*cos(coe_pl(6)*pi/180) ) ;

% orbital parameters of arrival trajectory
h1=coe_in(1); % specific angular momentum
e1=coe_in(2); % eccentricity
TA_in = coe_in(6)  ;

v_ort_1 = (muSun/h1) * (1 + e1 * cos(TA_in) ) ; % orthogonal component of spacecraft heliocentric velocity at arrival on planet's SOI
v_radial_1 = (muSun/h1) * e1 * sin(TA_in) ; % radial component of spacecraft heliocentric velocity at arrival on planet's SOI
% heliocentric spacecraft's velocity vector in vs frame(v=direction of
% heliocentric velocity of the planet and s=direction to sun)
v_in_vs = [ v_ort_1 , -v_radial_1 , 0];
% angle between V_in and v versor
alpha_in = atan2(v_in_vs(2),v_in_vs(1)) ; % equal to -path_angle
% excess velocity vector in vs frame at entrance in SOI
v_inf_in = v_in_vs - [ norm(V) , 0 ,0 ] ;
v_inf = norm(v_inf_in);

% orbital elements of flyby hyperbola
e = 1 + rp * v_inf^2 / mu ; % eccentricity
h = rp * sqrt( v_inf^2 + 2*mu / rp ) ; % specific angular momentum
TA_inf = acos(-1/e) ; % asyptotic value of the true anomaly
beta = pi - TA_inf ; % angle between asymptote and apse line
delta = 2*asin(1/e); % turn angle
a = mu*v_inf^2 ;  % semimajor axis
aiming = a * sqrt(e^2 -1 )  ; % aiming radius

TA_in = acos ( e^-1 * (h^2/(mu*SOI) -1 ) ) ; % true anomaly at SOI entrance
F_in = etheta2E(e, TA_in); % eccentric anomaly at SOI entrance
M_in = e * sinh(F_in) - F_in ; % mean anomaly at SOI entrance
n_in = mu^2 * (e^2-1) ^1.5 /h^3 ; % mean motion at SOI entrance
TOF_in = M_in / n_in ; % TOF from SOI entrance to periapsis
TOF_fb = 2*TOF_in ; % total flyby's TOF

% angle between v_inf_in and v axis
phi_in = atan2(v_inf_in(2),v_inf_in(1) ) ;

% Assuming a trailing side flyby,in order to improve heliocentric speed of
% the spacecraft :
if phi_in<0
    phi_out = phi_in + delta ; % angle between v_inf_out and v axis
    w= beta + phi_out +pi ; % argument of perigeo(from v axis to perifocal p axis)
else
    phi_out = phi_in - delta ;  % angle between v_inf_out and v axis
    w= phi_out + TA_inf ; % argument of perigeo(from v axis to perifocal p axis)
end

% angle between v versor and apsis line (argument of hyperbola's periapsis in vs frame)
RA = 0; 
i = 0; 
R_from_per =  (R_to_perifocal(RA,i,w))' ; % rotation matrix from hyperbola's perifocal frame to vs frame

r_periapsis_vs = R_from_per * [rp,0,0]'  ;
vp = h / rp ; % speed at hyperbola's periapsis
v_periapsis_vs = R_from_per * [0,vp,0]' ; % periapsis velocity in vs frame
 
% excess velocity vector in vs frame at exit from SOI
v_inf_out = [ v_inf*cos(phi_out) , v_inf*sin(phi_out) , 0 ];
% spacecraft heliocentric velocity vector at exit from SOI in vs frame
v_out_vs = v_inf_out + [norm(V),0,0] ;
v_ort_2 = v_out_vs(1); % exit velocity's orthogonal component
v_radial_2 = -v_out_vs(2); % exit velocity's radial component
 
% orbital elements of exit trajectory
h2 = norm(R) * v_ort_2 ; % specific angular momentum(km^2/s)

e2sinTh2 =  v_radial_2 * h2 / muSun ; 
e2cosTh2 =  (h2^2/muSun/norm(R) -1 ) ; 
TA_out = atan2(e2sinTh2 , e2cosTh2) ; % true anomaly in post flyby orbit 

e2 = e2cosTh2 / cos(TA_out) ; % post flyby orbit 's eccentricity
a2 = h2^2/muSun/(1-e2^2) ; % post flyby orbit 's semimajor axis
RA_2 = 0; 
incl_2 = 0; 
w_2 = (coe_pl(5) + coe_pl(6)- TA_out/deg ) *deg; % argument of perihelium(rad)

coe_out = [h2 e2 RA_2 incl_2 w_2 TA_out a2] ; % post flyby orbital elements

% angle between planet velocity vector and x axis
ang1 = atan2(V(2),V(1));
% angle between out spacecraft velocity vector and planet velocity vector
alpha_out = atan2(v_out_vs(2),v_out_vs(1)) ;
% angle between out spacecraft velocity vector and x axis
ang3 = ang1+alpha_out ;
v_out_xy = norm(v_out_vs) * [cos(ang3),sin(ang3),0]; % rotate v_out in xy frame

%dVgained= 2*v_inf/e ; 
dV = v_inf_out-v_inf_in;
dVgained = norm(dV) ; % deltaV gained with the gravity assist

if graphic==true
    % Create FlyBy plot
    figure
    hold all
    whitebg('k')
    
    if fb_number==1 % Venus flyby
    
        xlim([-1e5  1e5 ])
        ylim([-1e5  1e5 ])
        axis equal

        title(['Venus flyby (',num2str(date(3)),'/',num2str(date(2)),'/',num2str(date(1)),')'])

        scale_vel = 5e3 ; % scale factor to visualize velocity vectors 
 
        % plot flyby hyperbola
        OrbPlot(r_periapsis_vs,v_periapsis_vs,mu,'r',2,TA_in)
        OrbPlot(r_periapsis_vs,v_periapsis_vs,mu,'r',2,-TA_in)

        scatter3(0,0,0,100,'w','filled') % flyby planet

        quiver3(5e4,0,0,5e4,0,0,'color','b')  % uv versor(planet velocity direction)
        text(8e4,6e3,'u_{v}')

        quiver3(0,1e4,0,0,4e4,0,'color','y') % us versor(to the Sun)
        text(5e3,3e4,'u_{s}')

        % arrival relative velocity
        quiver3(8e4,8e4,0,v_inf_in(1)*scale_vel,v_inf_in(2)*scale_vel,v_inf_in(3)*scale_vel,'color','r','linewidth',1)
        text(9e4,7e4,'v_{inf} in')
        % Planet velocity
        quiver3(8e4-norm(V)*0.9*scale_vel,8e4,0,norm(V)*scale_vel,0,0,'color','b')
        text(0,9e4,'V_{Venus}')
        % Heliocentric arrival velocity of the spacecraft
        quiver3(8e4-norm(V)*0.9*scale_vel,8e4,0,v_in_vs(1)*scale_vel,v_in_vs(2)*scale_vel,0,'color','g')
        text(-4e4,6e4,['V_{in} (',num2str(norm(v_in_vs)),' Km/s )'])

        % exit relative velocity
        quiver3(8e4,-6e4,0,v_inf_out(1)*scale_vel,v_inf_out(2)*scale_vel,v_inf_out(3)*scale_vel,'color','r','linewidth',1)
        text(9e4,-6.5e4,'v_{inf} out')
        % Planet velocity
        quiver3(8e4-norm(V)*0.9*scale_vel,-6e4,0,norm(V)*scale_vel,0,0,'color','b')
        text(-5e4,-5e4,'V_{Venus}')
        % Heliocentric exit velocity of the spacecraft
        quiver3(8e4-norm(V)*0.9*scale_vel,-6e4,0,v_out_vs(1)*scale_vel,v_out_vs(2)*scale_vel,0,'color','g')
        text(-3e4,-8e4,['V_{out} (',num2str(norm(v_out_vs)),' Km/s )'])
    
        % triangle of relative velocity vectors showing the deltaV vector
        quiver3(-1e5,0,0,v_inf_in(1)*scale_vel,v_inf_in(2)*scale_vel,v_inf_in(3)*scale_vel,'color','r','linewidth',1)
        quiver3(-1e5,0,0,v_inf_out(1)*scale_vel,v_inf_out(2)*scale_vel,v_inf_out(3)*scale_vel,'color','r','linewidth',1)
        quiver3(-1e5+v_inf_in(1)*0.9*scale_vel,v_inf_in(2)*0.9*scale_vel,v_inf_in(3)*0.9*scale_vel,dV(1)*scale_vel,dV(2)*scale_vel,dV(3)*scale_vel,'color','w','linewidth',1)
        text(-8.5e4,-3e4,['dV (',num2str(dVgained),' Km/s )'])
        
        % turn angle
        ang([-1e5,0],1.5e4,[phi_in phi_out],'w') % turn angle
        ang([-1e5,0],1.52e4,[phi_in phi_out],'w') 
        ang([-1e5,0],1.54e4,[phi_in phi_out],'w')        
        text(-9e4,-0.5e4,['\delta =',num2str(delta/deg),'\circ'])

        
    elseif fb_number==2 % First Earth FlyBy
        
        xlim([-1e5  1e5 ])
        ylim([-1e5  1e5 ])
        axis equal

        title(['First Earth flyby (',num2str(date(3)),'/',num2str(date(2)),'/',num2str(date(1)),')'])

        scale_vel = 5e3 ; % scale factor to visualize velocity vectors 

        % plot flyby hyperbola
        OrbPlot(r_periapsis_vs,v_periapsis_vs,mu,'r',2,TA_in)
        OrbPlot(r_periapsis_vs,v_periapsis_vs,mu,'r',2,-TA_in)

        scatter3(0,0,0,100,'w','filled') % flyby planet

        quiver3(5e4,0,0,5e4,0,0,'color','b')  % uv versor(planet velocity direction)
        text(8e4,6e3,'u_{v}')

        quiver3(0,6e4,0,0,4e4,0,'color','y') % us versor(to the Sun)
        text(5e3,8e4,'u_{s}')

        % arrival relative velocity
        quiver3(10e4,-8e4,0,v_inf_in(1)*scale_vel,v_inf_in(2)*scale_vel,v_inf_in(3)*scale_vel,'color','r','linewidth',1)
        text(10e4,-6e4,'v_{inf} in')
        % Planet velocity
        quiver3(10e4-norm(V)*0.9*scale_vel,-8e4,0,norm(V)*scale_vel,0,0,'color','b')
        text(5e4,-9e4,'V_{Earth}')
        % Heliocentric arrival velocity of the spacecraft
        quiver3(10e4-norm(V)*0.9*scale_vel,-8e4,0,v_in_vs(1)*scale_vel,v_in_vs(2)*scale_vel,0,'color','g')
        text(1e4,-4e4,['V_{in} (',num2str(norm(v_in_vs)),' Km/s )'])

        % exit relative velocity
        quiver3(10e4,4e4,0,v_inf_out(1)*scale_vel,v_inf_out(2)*scale_vel,v_inf_out(3)*scale_vel,'color','r','linewidth',1)
        text(12e4,5e4,'v_{inf} out')
        % Planet velocity
        quiver3(10e4-norm(V)*0.9*scale_vel,4e4,0,norm(V)*scale_vel,0,0,'color','b')
        text(5e4,3e4,'V_{Earth}')
        % Heliocentric exit velocity of the spacecraft
        quiver3(10e4-norm(V)*0.9*scale_vel,4e4,0,v_out_vs(1)*scale_vel,v_out_vs(2)*scale_vel,0,'color','g')
        text(6e4,7.5e4,['V_{out} (',num2str(norm(v_out_vs)),' Km/s )'])
    
        % triangle of relative velocity vectors showing the deltaV vector
        quiver3(-5e4,-2e4,0,v_inf_in(1)*scale_vel,v_inf_in(2)*scale_vel,v_inf_in(3)*scale_vel,'color','r','linewidth',1)
        quiver3(-5e4,-2e4,0,v_inf_out(1)*scale_vel,v_inf_out(2)*scale_vel,v_inf_out(3)*scale_vel,'color','r','linewidth',1)
        quiver3(-5e4+v_inf_in(1)*0.9*scale_vel,-2e4+v_inf_in(2)*0.9*scale_vel,v_inf_in(3)*0.9*scale_vel,dV(1)*scale_vel,dV(2)*scale_vel,dV(3)*scale_vel,'color','w','linewidth',1)
        text(-5e4,3e4,['dV (',num2str(dVgained),' Km/s )'])
        
        % turn angle
        ang([-5e4,-2e4],1.5e4,[phi_in phi_out],'w') 
        ang([-5e4,-2e4],1.52e4,[phi_in phi_out],'w') 
        ang([-5e4,-2e4],1.54e4,[phi_in phi_out],'w')        
        text(-4e4,-2e4,['\delta =',num2str(delta/deg),'\circ'])
               
    elseif fb_number==3 % Second Earth Flyby
    
        xlim([-1e5  1e5 ])
        ylim([-1e5  1e5 ])
        %axis equal

        title(['Second Earth flyby (',num2str(date(3)),'/',num2str(date(2)),'/',num2str(date(1)),')'])

        scale_vel = 5e3 ; % scale factor to visualize velocity vectors 

        % plot flyby hyperbola
        OrbPlot(r_periapsis_vs,v_periapsis_vs,mu,'r',2,TA_in)
        OrbPlot(r_periapsis_vs,v_periapsis_vs,mu,'r',2,-TA_in)

        scatter3(0,0,0,100,'w','filled') % flyby planet

        quiver3(5e4,0,0,5e4,0,0,'color','b')  % uv versor(planet velocity direction)
        text(8e4,6e3,'u_{v}')

        quiver3(0,6e4,0,0,4e4,0,'color','y') % us versor(to the Sun)
        text(5e3,8e4,'u_{s}')

        % arrival relative velocity
        quiver3(5e4,-8e4,0,v_inf_in(1)*scale_vel,v_inf_in(2)*scale_vel,v_inf_in(3)*scale_vel,'color','r','linewidth',1)
        text(7e4,-7e4,'v_{inf} in')
        % Planet velocity
        quiver3(5e4-norm(V)*0.9*scale_vel,-8e4,0,norm(V)*scale_vel,0,0,'color','b')
        text(0,-9e4,'V_{Earth}')
        % Heliocentric arrival velocity of the spacecraft
        quiver3(5e4-norm(V)*0.9*scale_vel,-8e4,0,v_in_vs(1)*scale_vel,v_in_vs(2)*scale_vel,0,'color','g')
        text(-2e4,-5e4,['V_{in} (',num2str(norm(v_in_vs)),' Km/s )'])

        % exit relative velocity
        quiver3(5e4,4e4,0,v_inf_out(1)*scale_vel,v_inf_out(2)*scale_vel,v_inf_out(3)*scale_vel,'color','r','linewidth',1)
        text(6e4,3e4,'v_{inf} out')
        % Planet velocity
        quiver3(5e4-norm(V)*0.9*scale_vel,4e4,0,norm(V)*scale_vel,0,0,'color','b')
        text(0,3e4,'V_{Earth}')
        % Heliocentric exit velocity of the spacecraft
        quiver3(5e4-norm(V)*0.9*scale_vel,4e4,0,v_out_vs(1)*scale_vel,v_out_vs(2)*scale_vel,0,'color','g')
        text(0e4,5e4,['V_{out} (',num2str(norm(v_out_vs)),' Km/s )'])
    
        % triangle of relative velocity vectors showing the deltaV vector
        quiver3(-7e4,0,0,v_inf_in(1)*scale_vel,v_inf_in(2)*scale_vel,v_inf_in(3)*scale_vel,'color','r','linewidth',1)
        quiver3(-7e4,0,0,v_inf_out(1)*scale_vel,v_inf_out(2)*scale_vel,v_inf_out(3)*scale_vel,'color','r','linewidth',1)
        quiver3(-7e4+v_inf_in(1)*0.9*scale_vel,v_inf_in(2)*0.9*scale_vel,v_inf_in(3)*0.9*scale_vel,dV(1)*scale_vel,dV(2)*scale_vel,dV(3)*scale_vel,'color','w','linewidth',1)
        text(-3e4,2e4,['dV (',num2str(dVgained),' Km/s )'])
        
        % turn angle
        ang([-7e4,-0e4],1.5e4,[phi_in phi_out],'w') 
        ang([-7e4,-0e4],1.52e4,[phi_in phi_out],'w')
        ang([-7e4,-0e4],1.54e4,[phi_in phi_out],'w')
        text(-6e4,0e4,['\delta =',num2str(delta/deg),'\circ'])

    end
    
end