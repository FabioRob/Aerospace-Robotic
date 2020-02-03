function [date_next_target , zp_fb ] = fb_date_and_altitude ( fb_planet_id , target_planet_id , date_current_fb , date_init_target, date_max_target , zp_min , zp_inc , coe_in ) 

% -----------------------------------------------------------------------%
%
% This function computes altitude's flyby and date in which the post-flyby target planet will be reached
% If a solution is found, it outputs that, else, the output corresponds to
% the neaerest distance spacecraft-target planet obtained
%
% Arguments:
%
% fb_planet_id      - planet_id of flyby planet
% target_planet_id  - planet_id of target planet
% date_current_fb   - date of flyby
% date_init_target  - date from which a possible encounter date with the target planet must be found
% date_max_target   - date until which a possible encounter date with the target planet must be found
% zp_min            - minimum altitude of the flyby
% zp_inc            - percentual of SOI that define the increment for zp(e.g. 1e-4)
% coe_in            - orbital elements of arrival heliocentric orbit
%
% Outputs :
%
% date_next_target  - date in which encounter with next target will happen(or in which there is minimum distance)
% zp                - altitude of periapsis flyby hyperbola
%
% -----------------------------------------------------------------------%


muSun=132712439935; % sun's gravitational parameter
[~,radius,SOI] = physical_data(fb_planet_id);

[~, ra_fb_planet, ~ , ~ ] = planet_elements_and_sv(3,fb_planet_id, date_current_fb(1),date_current_fb(2),date_current_fb(3),date_current_fb(4), ...
    date_current_fb(5),date_current_fb(6));

jd_fb = date2JD(date_current_fb(1),date_current_fb(2),date_current_fb(3),date_current_fb(4),date_current_fb(5),date_current_fb(6));
jd_target = date2JD( date_init_target(1),date_init_target(2),date_init_target(3),date_init_target(4),date_init_target(5),date_init_target(6));

zp_increment = SOI * zp_inc  ;
zp = (zp_min:zp_increment:(SOI-radius)) ;

jd_target_max = date2JD( date_max_target(1),date_max_target(2),date_max_target(3),date_max_target(4),date_max_target(5),date_max_target(6));

jd = (jd_target:1:jd_target_max) ; % test dates in given range

min_distance = 1e12 ;

for i = 1:length(zp)
    for j = 1 : length(jd)
        
        [ ~ , ~ , v_out ] = FlyBy ( fb_planet_id , date_current_fb , coe_in  , zp(i) ,false,0 ); % flyby

        [year, month, day, hour, minute, second] = JD2date(jd(j) ) ;
        date_next_target = [year, month, day, hour, minute, second] ;
                
        tof = (jd(j) - jd_fb ) * 86400 ; % sec
        
        [r_s_arr_target , ~ ] = rv_from_r0v0(ra_fb_planet, v_out , tof , muSun) ;
        
        [~, ra_E_2, ~, ~ ] = planet_elements_and_sv(3,target_planet_id , date_next_target(1),date_next_target(2),date_next_target(3),date_next_target(4),date_next_target(5),date_next_target(6));
        distance_spacecraft_target = norm(r_s_arr_target - ra_E_2 ) ;  % it must be less then SOI
        
        if distance_spacecraft_target<SOI
            zp_fb = zp(i) ;
            disp('got solution')
            return
        end
        
        if distance_spacecraft_target<min_distance
            zp_fb = zp(i) ;
            min_distance = distance_spacecraft_target;
            date_min = date_next_target ;
        end
    end
end
disp('most near solution ... ')
date_next_target = date_min ;

    
    
    