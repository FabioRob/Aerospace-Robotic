% This script creates spacecraft's trajectory animation 

load 'state_vectors_for_animation'
load 'dates_for_animation'

% To create 4 separates animations :

%{
spacecraft_animation(3, 3 , 2 ,date_dep, date_fb_Venus, rd_s_E , vd_s_E , 'E2V' )
spacecraft_animation(3, 2 , 3 ,date_fb_Venus, date_fb_Earth_1, ra_s_V ,v_out_1 , 'V2E')
spacecraft_animation(3, 3 , 3 ,date_fb_Earth_1, date_fb_Earth_2, ra_E_1 ,v_out_2 , 'E2E')
spacecraft_animation(3, 3 , 7 ,date_fb_Earth_2, date_Uranus, ra_E_2 ,vd_s , 'E2U')
%}

% To create one animation with the entire spacecraft's trajectory

animation_dir = 'Trajectory' ;

muSun=132712439935; % sun's gravitational parameter
    
JD_dep = date2JD(date_dep(1),date_dep(2),date_dep(3),date_dep(4),date_dep(5),date_dep(6));
JD_arr_V = date2JD(date_fb_Venus(1),date_fb_Venus(2),date_fb_Venus(3),date_fb_Venus(4),date_fb_Venus(5),date_fb_Venus(6));
JD_arr_E_1 = date2JD(date_fb_Earth_1(1),date_fb_Earth_1(2),date_fb_Earth_1(3),date_fb_Earth_1(4),date_fb_Earth_1(5),date_fb_Earth_1(6));
JD_arr_E_2 = date2JD(date_fb_Earth_2(1),date_fb_Earth_2(2),date_fb_Earth_2(3),date_fb_Earth_2(4),date_fb_Earth_2(5),date_fb_Earth_2(6));
JD_arr_U = date2JD(date_Uranus(1),date_Uranus(2),date_Uranus(3),date_Uranus(4),date_Uranus(5),date_Uranus(6));

tof_E_V = (JD_arr_V-JD_dep)*86400;
tof_V_E = (JD_arr_E_1-JD_arr_V)*86400;
tof_E_E = (JD_arr_E_2-JD_arr_E_1)*86400;
tof_E_U = (JD_arr_U - JD_arr_E_2)*86400;


[rf_V,vf_V] = rv_from_r0v0( rd_s_E ,vd_s_E,tof_E_V,muSun); % spacecraft state vector on Venus arrival
coe_arr_V = coe_from_sv(rf_V,vf_V,muSun);
TA_arr_V = coe_arr_V(6);

[rf_E_1,vf_E_1] = rv_from_r0v0( ra_s_V ,v_out_1 ,tof_V_E ,muSun); % spacecraft state vector on first Earth arrival
coe_arr_E_1 = coe_from_sv(rf_E_1,vf_E_1,muSun);
TA_arr_E_1 = coe_arr_E_1(6);

[rf_E_2,vf_E_2] = rv_from_r0v0( ra_E_1 , v_out_2 ,tof_E_E ,muSun); % spacecraft state vector on second Earth arrival
coe_arr_E_2 = coe_from_sv(rf_E_2,vf_E_2,muSun);
TA_arr_E_2 = coe_arr_E_2(6);

[rf_U,vf_U] = rv_from_r0v0( ra_E_2 , vd_s  ,tof_E_U ,muSun); % spacecraft state vector on Uranus arrival
coe_arr_U = coe_from_sv(rf_U,vf_U,muSun);
TA_arr_U = coe_arr_U(6);

frame_rate = 25 ;
mkdir(animation_dir);
workingDir=animation_dir;
h=figure('visible', 'off');
whitebg('k')
k = 0;
outputVideo = VideoWriter(fullfile(workingDir,[animation_dir,'movie.avi']));
outputVideo.FrameRate = frame_rate;
open(outputVideo)
hold all
view([0,90])

draw_sphere([0,0,0], 30, 'y', 1);

plotPlanetOrbit (3 , 2 , date_dep , 'g' , 1)
plotPlanetOrbit (3 , 3 , date_dep , 'w' , 1)
plotPlanetOrbit (3 , 4 , date_dep , 'r' , 1)
plotPlanetOrbit (3 , 5 , date_dep , 'b' , 1)
plotPlanetOrbit (3 , 6 , date_dep , 'y' , 1)
plotPlanetOrbit (3 , 7 , date_dep , 'c' , 1)

orbcol = [1,0.5,0] ;
OrbPlot(rd_s_E ,vd_s_E,muSun,orbcol,1,TA_arr_V)
OrbPlot(ra_s_V ,v_out_1,muSun,orbcol,1,TA_arr_E_1)
OrbPlot(ra_E_1 , v_out_2,muSun,orbcol,1,TA_arr_E_2)
OrbPlot(ra_E_2 , vd_s ,muSun,orbcol,1,TA_arr_U)

xlim([-1e9,1e9])
ylim([-0.5e9,3e9])
axis equal


[k1] = create_frames (outputVideo , h , frame_rate, JD_dep ,     JD_arr_V ,  rd_s_E ,vd_s_E ,  workingDir , 0  )  ; % E2V
[k2] = create_frames (outputVideo , h , frame_rate, JD_arr_V ,   JD_arr_E_1 ,ra_s_V ,v_out_1 , workingDir , k1  ) ; % V2E
[k3] = create_frames (outputVideo , h , frame_rate, JD_arr_E_1 , JD_arr_E_2 ,ra_E_1 , v_out_2, workingDir , k2  ) ; % E2E
[k4] = create_frames (outputVideo , h , frame_rate, JD_arr_E_2 , JD_arr_U ,  ra_E_2 , vd_s  ,  workingDir , k3  ) ; % E2U
close(outputVideo)



