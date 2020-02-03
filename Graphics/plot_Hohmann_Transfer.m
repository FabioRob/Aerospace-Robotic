function [V_E_H , V_U_H ] =  plot_Hohmann_Transfer(date_dep,date_arr)

% -----------------------------------------------------------------------%
%
% plot_Hohmann_Transfer function plots the Hohmann transfer ellipse from Earth to Uranus
% and computes departure and arrival velocity vectors
%
% Arguments :
%
% date_dep - Earth departure date
% date_arr - Uranus arrival date
%
% Outputs :
%
% V_E_H    - Spacecraft's velocity vector in heliocentric frame at departure from Earth
% V_U_H    - Spacecraft's velocity vector in heliocentric frame at arrival on Uranus
%
% -----------------------------------------------------------------------%

deg = pi/180; % to convert from degrees to radians
muSun= 132712439935; % Sun's gravitational parameter(km^3/s^2)

[ ~ , v_dep , ~ , v_arr , ~ , ~ , ~ , TOF_E_U  ] = Compute_Direct_Hohmann ( 3 , 7  );

figure
hold all
whitebg('k')    
title('Direct Hohmann Transfer from Earth to Uranus')

scatter(0,0,100,'y','filled') % Sun

xlim([-3.3e9 , 4e9])
ylim([-2.8e9,2.8e9])
grid on
axis equal

plotPlanetOrbit (3 , 4 , date_dep , 'r' , 1 ) % Mars orbit
plotPlanetOrbit (3 , 5 , date_dep , 'b' , 1 ) % Jupiter orbit
plotPlanetOrbit (3 , 6 , date_dep , 'y' , 1 ) % Saturn orbit

vel_scale_H = 1e8 ; % scale factor used to plot velocity vectors 

% Earth at departure
[oe_E_dep, R_E_dep, V_E_dep, ~] = ComputeAndPlotPlanetPosition (2,3,date_dep,'o','w');
text(0.04e9,-0.2e9,'\uparrow')
text(-0.2e9 , -0.3e9,[ 'Earth at departure date (',num2str(date_dep(3)),'/',num2str(date_dep(2)),'/',num2str(date_dep(1)),')'])

OrbPlot(R_E_dep,V_E_dep,muSun,'w',2)  % Earth orbit
text(-9e8,-1e8,' Earth orbit\rightarrow')

[~, ru, vu, ~] = planet_elements_and_sv(2,7, date_dep(1) ,date_dep(2), date_dep(3), date_dep(4), date_dep(5), date_dep(6)) ;
OrbPlot(ru, vu,muSun,'c',2)  % Uranus orbit
text(-2.5e9,-2.5e9,'Uranus orbit\rightarrow')

% Compute rotation matrix from perifocal hohmann transfer frame to
% heliocentric frame
RA_h = 0; % right ascension of transfer ellipse is zero
i_h = 0; % inclination of transfer ellipse is zero
w_h = oe_E_dep(6)*deg ; % compute argument of perihelium through Earth position vector at departure( angle from X axis to perifocal x axis)
R_from_per_hohmann =  (R_to_perifocal(RA_h,i_h,w_h))' ; % rotation matrix from perifocal frame to heliocentric frame

% velocity vectors at departure:
% Earth veocity vector
quiver3(R_E_dep(1),R_E_dep(2),R_E_dep(3),V_E_dep(1)*vel_scale_H,V_E_dep(2)*vel_scale_H,V_E_dep(3)*vel_scale_H,'Color','w','LineWidth',1)
text(0.5e9,0.65e9 , ['(',num2str(norm(V_E_dep),2),' km/s)'] )
text(0.5e9,0.5e9,'Earth velocity  \rightarrow')
% Spacecraft velocity vector in the Hohmann transfer trajectory at departure from Earth
% in heliocentric frame
V_E_H =( R_from_per_hohmann * [0,v_dep,0]' )';
quiver3(R_E_dep(1)+1e8,R_E_dep(2)-1e8*cos(w_h),R_E_dep(3),V_E_H(1)*vel_scale_H,V_E_H(2)*vel_scale_H,0,'Color','g','LineWidth',1)
text(3e9,0.7e9,'\uparrow')
text(2.4e9,0.5e9,'velocity at departure')
text(2.4e9,0.4e9,'in the Hohmann transfer')
text(2.6e9,0.3e9,['(',num2str(norm(V_E_H),2),' km/s)'] )

% Uranus at departure from Earth
ComputeAndPlotPlanetPosition (2,7,date_dep,'o','w');
text(2.12e9,1.8e9,'\uparrow')
text(1.2e9,1.6e9,['Uranus at departure date( ',num2str(date_dep(3)),'/',num2str(date_dep(2)),'/',num2str(date_dep(1)),')'])

% Plot planets positions at arrival date
[~, ~, ~, ~ ]=ComputeAndPlotPlanetPosition(2,3,date_arr,'*','w');
text(0.1e9,-0.1e9,['\leftarrow Earth at arrival date( ',num2str(date_arr(3)),'/',num2str(date_arr(2)),'/',num2str(date_arr(1)),')'])

[~, R_U_arr, V_U_arr,~]=ComputeAndPlotPlanetPosition(2,7,date_arr,'*','w');
text(-0.8e9,2.7e9,['\leftarrow Uranus at arrival date( ',num2str(date_arr(3)),'/',num2str(date_arr(2)),'/',num2str(date_arr(1)),')'])

% velocity vectors at arrival:
% Uranus velocity vector
quiver3(R_U_arr(1),R_U_arr(2),R_U_arr(3),V_U_arr(1)*vel_scale_H,V_U_arr(2)*vel_scale_H,V_U_arr(3)*vel_scale_H,'Color','c','LineWidth',2)
text(-1.5e9,2.4e9,'\uparrow')
text(-2e9,2.2e9,'Uranus velocity')
text(-2e9, 2e9, ['(' ,num2str(norm(V_U_arr),2),' km/s)'])
% Spacecraft velocity vector in the Hohmann transfer trajectory at arrival
% on Uranus, in heliocentric frame
V_U_H =( R_from_per_hohmann * [0, -v_arr , 0]' )' ;
quiver3(R_U_arr(1)-1e8,R_U_arr(2)+1e8*cos(w_h),R_U_arr(3),V_U_H(1)*vel_scale_H,V_U_H(2)*vel_scale_H,0,'Color','g','LineWidth',1)
text(-2.5e9,2.8e9,'velocity at arrival in the')
text(-2.5e9,2.7e9,'Hohmann transfer \rightarrow')
text(-2.5e9,2.6e9,['(' ,num2str(v_arr,2),' km/s)'])


% Plot Hohmann transfer ellipse apse line
Hohmann_apseLine = [R_E_dep ; R_U_arr ];
plot3(Hohmann_apseLine(:,1),Hohmann_apseLine(:,2),Hohmann_apseLine(:,3),'--')

OrbPlot(R_E_dep, V_E_H  ,muSun, 'g' , 2 , pi) % plot Hohmann transfer trajectory
text(0,2e9,'Hohmann transfer')
text(0.1e9,1.9e9,['(TOF\sim', num2str(TOF_E_U/86400/365,2),' years)'])
