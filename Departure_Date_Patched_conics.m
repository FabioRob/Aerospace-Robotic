% This script creates "pork chops" for transfer Earth-Venus
% figure 1(C3L-v_inf_arr) : contour plot with caracteristic energy(red) 
% C3L(Km^2/s^2 , equal to the square of spacecraft's relative velocity at departure) 
% and spacecraft's relative velocity on arrival(blue) v_inf_arr(Km/s)
% figure 2(C3L-TOF) : C3L(red) and time of flight(blu,days)

load_data

row = 0;
col = 0;

fprintf('\nprogram porkchop\n');
fprintf('\n< interplanetary pork chop plots >');

month= 1 ;
day= 1 ;
year = 2020; 
fprintf('\n\n\nnominal launch date\n');
start_date_E =[ month,day,year ]
% launch must occurs after 1st january 2020(progect specification 2)
JD_dep = julian(month, day, year);

month= 4 ;
day= 1 ;
year = 2020; 
fprintf('\n\nnominal arrival date\n');
start_date_V =[ month,day,year ]
JD_arr = julian(month, day, year);


fprintf('\n\nlaunch date span in days\n');
spanx = 365 


fprintf('\n\narrival date span in days\n');
spany = 365 
  


fprintf('\n\nstep size in days\n');
step = 5
   
c3_levels = [5,8,9,10,12,15,16,20,30,40,50];   

vinf_levels = [4.0,4.5,5.0,5.3,5.6,6.0,6.5,7.0];

tof_levels = [100,125,150,175,200,225,250,275,300,350,400,450];

dvt_levels = [5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12.0, 13.0, 14.0, 15.0];

fprintf('\n\nplease wait, computing data ....\n\n');

% compute C3L , v-infinity and TOF 
% for i = JD_dep - spanx : step: JD_dep + spanx
for i = JD_dep : step: JD_dep + spanx
    row = row + 1;
    col = 0;
    
    JDi = i;
    
    %for j = JD_arr - spany: step: JD_arr + spany
     for j = JD_arr : step: JD_arr + spany   
        col = col + 1;
        
        JDf = j;
        
        con_x(row, col) = JDi - JD_dep;
        con_y(row, col) = JDf - JD_arr;
                
        [dv1, dv2] = pcfunc (JDi, JDf);
                      
        c3l = norm(dv1) * norm(dv1);
        c3_dep(row, col) = c3l;
        vinf_arr(row, col) = norm(dv2);
               
        % compute flight time in days
        tof(row, col) = JDf - JDi;
                
    end
    
end

% create string representations of launch and arrival dates
[cdstr1, utstr1] = jd2str(JD_dep);
[cdstr2, utstr2] = jd2str(JD_arr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create combined C3l and arrival v-infinity plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
whitebg('w')

[c1, h1] = contour(con_x, con_y, c3_dep, c3_levels, 'r');

clabel(c1, h1);

hold on;

[c2, h2] = contour(con_x, con_y, vinf_arr, vinf_levels, 'b');
 
clabel(c2, h2);

title('Earth-to-Venus - C3L and arrival v-infinity', 'FontSize', 16');

xlabel(['days relative to Earth departure date of ', cdstr1]);
xlim([0,170])

ylabel(['days relative to Venus arrival date of ', cdstr2]);

legend([h1(1), h2(1)], 'C3L (km^2/sec^2)', 'arrival v_\infty (km/sec)', 'Location','NorthWest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create combined C3L and flight time plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
whitebg('w')

[c1, h1] = contour(con_x, con_y, tof, tof_levels, 'b');

clabel(c1, h1);

hold on;

[c2, h2] = contour(con_x, con_y, c3_dep, c3_levels, 'r');
 
clabel(c2, h2);

title('Earth-to-Venus - C3L and flight time', 'FontSize', 16);

xlabel(['days relative to Earth departure date of ', cdstr1]);

ylabel(['days relative to Venus arrival date of ', cdstr2]);

legend([h1(1), h2(1)], 'TOF (days)', 'C3L (km^2/sec^2)', 'Location','NorthWest');
