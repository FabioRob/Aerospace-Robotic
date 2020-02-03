function [ date_min_dist , zp_fb_min_dist , min_distance , date_min_burn , zp_fb_min_burn , min_burn ] = resolve_fb( fb_planet_id , target_planet_id , date_current_fb , date_init_target, date_max_target , zp_inc , coe_in ) 

% -----------------------------------------------------------------------%
%
% This function search best flyby altitudes(at a given flyby date
% date_current_fb) on fb_planet_id  and best dates to get to next target
% planet(after flyby) which minimize, respectively:
% the neaerest distance spacecraft-target planet obtained : date_next_target , zp_fb
% the additional burn required after flyby to put in the orbit to target planet: date_min_burn , zp_fb_min_burn
%
% Arguments:
%
% fb_planet_id      - planet_id of flyby planet
% target_planet_id  - planet_id of target planet
% date_current_fb   - date of flyby
% date_init_target  - date from which a possible encounter date with the target planet must be found
% date_max_target   - date until which a possible encounter date with the target planet must be found
% zp_inc            - percentual of SOI that define the increment for zp(e.g. 1e-4)
% coe_in            - orbital elements of arrival heliocentric orbit
%
% Outputs :
%
% date_min_dist     - date in which encounter with next target will happen(or in which there is minimum distance)
% zp_fb_min_dist    - altitude of periapsis flyby hyperbola in the above case
% min_distance      - minimum distance target-spacecraft
% date_min_burn     - date in which the additional burn required is minimum
% zp_fb_min_burn    - altitude of periapsis flyby in the above case
% min_burn          - minimum burn required to get to the target planet
%
% -----------------------------------------------------------------------%

muSun=132712439935; % sun's gravitational parameter
[~,radius,SOI] = physical_data(fb_planet_id);

[~, ra_fb_planet, ~ , ~ ] = planet_elements_and_sv(3,fb_planet_id, date_current_fb(1),date_current_fb(2),date_current_fb(3),date_current_fb(4), ...
    date_current_fb(5),date_current_fb(6));

jd_fb = date2JD(date_current_fb(1),date_current_fb(2),date_current_fb(3),date_current_fb(4),date_current_fb(5),date_current_fb(6));
jd_target = date2JD( date_init_target(1),date_init_target(2),date_init_target(3),date_init_target(4),date_init_target(5),date_init_target(6));

zp_increment = SOI * zp_inc  ;
zp = (0:zp_increment:(SOI-radius)) ;

jd_target_max = date2JD( date_max_target(1),date_max_target(2),date_max_target(3),date_max_target(4),date_max_target(5),date_max_target(6));

jd = (jd_target:1:jd_target_max) ; % test dates in given range

rd_s = zeros(length(jd),3) ;
vd_s = zeros(length(jd),3) ;

for k = 1:length(jd)
    [year, month, day, hour, minute, second] = JD2date(jd(k) ) ;
     date_arr = [year, month, day, hour, minute, second] ;
    [~ , rd ,vd , ~ ,~ ,~ , ~ ,~ ,~  ] = lambert_trajectory ( fb_planet_id , target_planet_id , date_current_fb , date_arr ,  false , 2) ;
    rd_s(k,:)=rd;
    vd_s(k,:)=vd;
end

coe_out=zeros(length(zp),7);
v_out=zeros(length(zp),3);

for r = 1:length(zp)
    [ ~ , coe , v ] = FlyBy ( fb_planet_id , date_current_fb , coe_in  , zp(r) ,false,0 ); % flyby
    coe_out(r,:) = coe;
    v_out(r,:) = v;
end

min_distance = 1e12 ;

min_burn = 1e10 ;


for r = 1:length(zp)
    for k = 1 : length(jd)                

        [year, month, day, hour, minute, second] = JD2date(jd(k) ) ;
        date_next_target = [year, month, day, hour, minute, second] ;
                
        tof = (jd(k) - jd_fb ) * 86400 ; % sec
        
        [r_s_arr_target , ~ ] = rv_from_r0v0(ra_fb_planet, v_out(r,:) , tof , muSun) ;
        
        [~, ra_pla_target, ~, ~ ] = planet_elements_and_sv(3,target_planet_id , date_next_target(1),date_next_target(2),date_next_target(3),date_next_target(4),date_next_target(5),date_next_target(6));
        distance_spacecraft_target = norm(r_s_arr_target - ra_pla_target ) ;   % it must be less then SOI
         
        coe_from_d2t = coe_from_sv(rd_s(k,:),vd_s(k,:),muSun); 

        path_angle_1 = atan2 ( coe_out(r,2)*sin(coe_out(r,6)) , 1 + coe_out(r,2)*cos(coe_out(r,6)) ) ; 
        path_angle_2 = atan2 ( coe_from_d2t(2)*sin(coe_from_d2t(6)) , 1 + coe_from_d2t(2)*cos(coe_from_d2t(6)) ) ; 

        dV_burn = sqrt(norm(v_out(r,:))^2 + norm(vd_s(k,:))^2 -2*norm(v_out(r,:))* norm(vd_s(k,:))*cos(path_angle_2-path_angle_1)) ;
                
        if dV_burn<min_burn
         
            date_min_burn=date_next_target;
            zp_fb_min_burn=zp(r);
            
            min_burn=dV_burn  ;    
        end
        
        if distance_spacecraft_target<min_distance
            
            zp_fb_min_dist = zp(r) ;
            date_min_dist=date_next_target;
            
            min_distance = distance_spacecraft_target;
        end
      
        
    end
    
end

