function [ dV_burn ] = deltaV_target (fb_planet_id ,date_fb , coe_out , v_out , target_planet_id , date_arr ) 

% This function computes deltaV burn required to change from post orbit
% flyby(on fb_planet_id , in date_fb , with post-flyby orbital elements coe_out and exit velocity v_out)
% to an orbit directed to the target planet on date_arr (found solving a lambert problem)

muSun=132712439935; % sun's gravitational parameter

% se volessi con un flyby sul pianeta fb_planet_id in date_fb, volessi direttamente arrivare su target_planet_id in data date_arr
[~ , rd_s ,vd_s , ~ ,~ ,~ , ~ ,~ ,~  ] = lambert_trajectory ( fb_planet_id , target_planet_id , date_fb , date_arr ,  false , 2) ;
coe_from_d2t = coe_from_sv(rd_s,vd_s,muSun); % elementi orbitali dell'orbita per il target planet , alla partenza dal pianeta di flyby

path_angle_1 = atan2 ( coe_out(2)*sin(coe_out(6)) , 1 + coe_out(2)*cos(coe_out(6)) ) ;
path_angle_2 = atan2 ( coe_from_d2t(2)*sin(coe_from_d2t(6)) , 1 + coe_from_d2t(2)*cos(coe_from_d2t(6)) ) ;

% quindi il deltaV necessario(pagina 325 Curtis) è :

dV_burn = sqrt(norm(v_out)^2 + norm(vd_s)^2 -2*norm(v_out)* norm(vd_s)*cos(path_angle_2-path_angle_1)) ; 
