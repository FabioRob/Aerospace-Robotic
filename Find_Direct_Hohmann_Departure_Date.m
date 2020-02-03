function [ date_dep , date_arr ] = Find_Direct_Hohmann_Departure_Date ( planet_id_dep , planet_id_arr , initial_date )

% -----------------------------------------------------------------------%
%
% Find_Direct_Hohmann_Departure_Date function iteratevely search for a date(starting from initial date)
% in which phase angle between the two planets is such that the spacecraft correctly rendezvous on arrival
% planet with an Hohmann transfer.
%
% Arguments :
%
% planet_id_dep - departure planet,from 1 to 7 (from Mercury to Uranus)
% planet_id_arr - arrival planet, from 1 to 7 (from Mercury to Uranus)
% initial_date  - launch date must be searched starting from this date
%                 ([year,month,day,hour,minute,second])
%
% Output :
%
% date_dep      - first possible date of departure for a Hohmann transfer
% date_arr      - arrival date
%
% -----------------------------------------------------------------------%


deg = pi/180; % to convert from degrees to radians
[ ~ , ~ , ~ , ~ , ~ , ~ , ~ , TOF_H  ] = Compute_Direct_Hohmann (planet_id_dep , planet_id_arr  ) ;

% Planets' orbit sideral periods(days)
orb_periods_planets = [ 87.97 , 224.7 , 365.256 , 1.881*365 , 11.86*365 , 29.46*365 , 84.01*365 ] ;

% angular velocity of arrival planet(rad/day)
n_arr = 2*pi / orb_periods_planets(planet_id_arr) ;

% Angular distance travelled by arrival planet(rad)
deltaTA_pl_arr = n_arr * TOF_H/86400;

% Initial phase between the two planets for a Hohmann transfer (rad)
phase_in = pi - deltaTA_pl_arr;

year=initial_date(1);
month=initial_date(2);
day=initial_date(3);
hour=initial_date(4);
minute=initial_date(5);
second=initial_date(6);
jd = date2JD(year, month, day, hour, minute, second);

phase = 0.0;
while ( abs(phase-phase_in) > 0.01 )
    
    [year, month, day, hour, minute, second] = JD2date(jd) ;  
    [oe_pl_dep_dep, ~, ~, ~ ] = planet_elements_and_sv(2,planet_id_dep, year, month, day, hour, minute, second);
    [oe_pl_arr_dep, ~ , ~,  ~] = planet_elements_and_sv(2,planet_id_arr, year, month, day, hour, minute, second);
    phase =zero_to_360(oe_pl_arr_dep(6) - oe_pl_dep_dep(6) ) *deg; % phase angle in radians in range [0-2*pi)
    
    jd = jd + 1 ;

end

jd_dep = jd-1;

% departure date
date_dep = [year, month, day, hour, minute, second] ; % first possible date for a Hohmann transfer

jd_arr = jd_dep + TOF_H/86400; % Julian day number of arrival date

% Arrival date
[year_arr, month_arr, day_arr, hour_arr, minute_arr, second_arr] = JD2date(jd_arr);
date_arr = [year_arr, month_arr, day_arr, hour_arr, minute_arr, second_arr] ;

