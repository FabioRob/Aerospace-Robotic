% interplanetary pork chop plots

% uses DE421 ephemeris

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all;

load_data

row = 0;
col = 0;

%clc; home;

fprintf('\nprogram porkchop\n');
fprintf('\n< interplanetary pork chop plots >');

fprintf('\n\n\nnominal launch date\n');
[month, day, year] = getdate;
JD_dep = julian(month, day, year);

fprintf('\n\nnominal arrival date\n');
[month, day, year] = getdate;
JD_arr = julian(month, day, year);

while(1)
   fprintf('\n\nplease input the launch date span in days\n');
   spanx = input('? ');
   
   if (spanx > 0)
      break;
   end
end 

while(1)
   fprintf('\n\nplease input the arrival date span in days\n');
   spany = input('? ');
   
   if (spany > 0)
      break;
   end
end

while(1)
   fprintf('\n\nplease input the step size in days\n');
   step = input('? ');
   
   if (step > 0)
      break;
   end
end

fprintf('\n\nplease input the launch energy contour levels in km^2/sec^2 ([] for defaults)\n');
c3_levels = input('? ');

if (size(c3_levels, 2) == 0)
   % c3_levels = [6,7,8,9,10,11,12,13,14,15,16,20,25,30,35,40,45,50];
   c3_levels = [5,8,9,10,12,15,16,20,30,40,50];
end
   
fprintf('\n\nplease input the arrival v-infinity contour levels in km/sec ([] for defaults)\n');
vinf_levels = input('? ');

if (size(vinf_levels, 2) == 0)
   vinf_levels = [4.0,4.2,4.5,4.8,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.5,7.0];
end

% fprintf('\n\nplease input the launch and arrival declination contour levels in degrees ([] for defaults)\n');
% dla_levels = input('? ');

% if (size(dla_levels, 2) == 0)
%    dla_levels = [-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30];
% end

% fprintf('\n\nplease input the launch and arrival right ascension contour levels in degrees ([] for defaults)\n');
% rla_levels = input('? ');

% if (size(rla_levels, 2) == 0)
%    rla_levels = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240, 255, 270, 285, 300, 315, 330];
% end

fprintf('\n\nplease input the time-of-flight contour levels in days ([] for defaults)\n');
tof_levels = input('? ');

if (size(tof_levels, 2) == 0)
   % tof_levels = [100, 150, 200, 250, 300, 350, 400, 450];
   tof_levels = [100,125,150,175,200,225,250,275,300,350,400,450];
end

fprintf('\n\nplease input the total delta-v contour levels in kilometers/second ([] for defaults)\n');
dvt_levels = input('? ');

if (size(dvt_levels, 2) == 0)
   dvt_levels = [5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12.0, 13.0, 14.0, 15.0];
end

fprintf('\n\nplease wait, computing data ....\n\n');

% compute C3L, v-infinity, DLA and RLA data
%for i = JD_dep - spanx : step: JD_dep + spanx
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
        
        % if (c3l <= 50)
        %    % compute orientation of the departure hyperbola
        %    decl_dep(row, col) = 90.0d0 - rtd * acos(dv1(3) / norm(dv1));
        %    rasc_dep(row, col) = rtd * atan3(dv1(2), dv1(1));
        % else
        %    decl_dep(row, col) = NaN;
        %    rasc_dep(row, col) = NaN;
        % end
       
        % compute arrival dla and rla
        % tmatrix = mme2000 (JD_dep + i);
        
        % dvmmee = tmatrix * (-dv2');
        
        % decl_arr(row, col) = 90.0d0 - rtd * acos(dvmmee(3) / norm(dv2));
        % rasc_arr(row, col) = rtd * atan3(dvmmee(2), dvmmee(1));
        
        % compute flight time in days
        tof(row, col) = JDf - JDi;
        
        % compute total delta-v in kilometers/second
        dvtotal(row, col) = norm(dv1) + norm(dv2);
        
    end
    
end

% create string representations of launch and arrival dates
[cdstr1, utstr1] = jd2str(JD_dep);
[cdstr2, utstr2] = jd2str(JD_arr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create combined C3l and arrival v-infinity plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);

[c1, h1] = contour(con_x, con_y, c3_dep, c3_levels, 'r');

clabel(c1, h1);

hold on;

[c2, h2] = contour(con_x, con_y, vinf_arr, vinf_levels, 'b');
 
clabel(c2, h2);

title('Earth-to-Venus - C3L and arrival v-infinity', 'FontSize', 16');

xlabel(['days relative to Earth departure date of ', cdstr1]);

ylabel(['days relative to Venus arrival date of ', cdstr2]);

legend([h1(1), h2(1)], 'C3L (km^2/sec^2)', 'arrival v_\infty (km/sec)', 'Location','NorthWest');

% print -depsc -tiff -r300 c3l_vinf.eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create combined C3L and DLA plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(2);

% [c1, h1] = contour(con_x, con_y, decl_dep, dla_levels, 'b');

% clabel(c1, h1);

% hold on;

% [c2, h2] = contour(con_x, con_y, c3_dep, c3_levels, 'r');
 
% clabel(c2, h2);

% title('Earth-to-Venus - Launch C3L and DLA', 'FontSize', 16');

% xlabel(['days relative to Earth departure date of ', cdstr1]);

% ylabel(['days relative to Venus arrival date of ', cdstr2]);

% legend([h1(1), h2(1)], 'DLA (deg)', 'C3L (km^2/sec^2)', 'Location','NorthWest');

% print -depsc -tiff -r300 c3l_dla.eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create combined C3L and RLA plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(3);

% [c1, h1] = contour(con_x, con_y, rasc_dep, rla_levels, 'b');

% clabel(c1, h1);

% hold on;

% [c2, h2] = contour(con_x, con_y, c3_dep, c3_levels, 'r');
 
% clabel(c2, h2);

% title('Earth-to-Venus - Launch C3L and RLA', 'FontSize', 16');

% xlabel(['days relative to Earth departure date of ', cdstr1]);

% ylabel(['days relative to Venus arrival date of ', cdstr2]);

% legend([h1(1), h2(1)], 'RLA (deg)', 'C3L (km^2/sec^2)', 'Location','NorthWest');

% print -depsc -tiff -r300 c3l_rla.eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create combined C3L and flight time plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);

[c1, h1] = contour(con_x, con_y, tof, tof_levels, 'b');

clabel(c1, h1);

hold on;

[c2, h2] = contour(con_x, con_y, c3_dep, c3_levels, 'r');
 
clabel(c2, h2);

title('Earth-to-Venus - C3L and flight time', 'FontSize', 16);

xlabel(['days relative to Earth departure date of ', cdstr1]);

ylabel(['days relative to Venus arrival date of ', cdstr2]);

legend([h1(1), h2(1)], 'TOF (days)', 'C3L (km^2/sec^2)', 'Location','NorthWest');

% print -depsc -tiff -r300 c3l_tof.eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create combined arrival v-infinity and DLA plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(5);

% [c1, h1] = contour(con_x, con_y, vinf_arr, vinf_levels, 'b');

% clabel(c1, h1);

% hold on;

% [c2, h2] = contour(con_x, con_y, decl_arr, dla_levels, 'r');
 
% clabel(c2, h2);

% title('Earth-to-Venus - arrival v-infinity and DLA', 'FontSize', 16);

% xlabel(['days relative to Earth departure date of ', cdstr1]);

% ylabel(['days relative to Venus arrival date of ', cdstr2]);

% legend([h1(1), h2(1)], 'v-infinity (km/sec)', 'DLA (degrees)', 'Location','NorthWest');

% print -depsc -tiff -r300 vinf_dla.eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create arrival v-infinity and RLA plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(6);

% [c1, h1] = contour(con_x, con_y, vinf_arr, vinf_levels, 'b');

% clabel(c1, h1);

% hold on;

% [c2, h2] = contour(con_x, con_y, rasc_arr, rla_levels, 'r');
 
% clabel(c2, h2);

% title('Earth-to-Venus - arrival v-infinity and RLA', 'FontSize', 16);

% xlabel(['days relative to Earth departure date of ', cdstr1]);

% ylabel(['days relative to Venus arrival date of ', cdstr2]);

% legend([h1(1), h2(1)], 'v-infinity (km/sec)', 'RLA (degrees)', 'Location','NorthWest');

% print -depsc -tiff -r300 vinf_rla.eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create toal delta-v and flight time plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure(7);

%[c1, h1] = contour(con_x, con_y, tof, tof_levels, 'b');

%clabel(c1, h1);

%hold on;

%[c2, h2] = contour(con_x, con_y, dvtotal, dvt_levels, 'r');
 
%clabel(c2, h2);

%title('Earth-to-Venus - total delta-v and flight time', 'FontSize', 16);

%xlabel(['days relative to Earth departure date of ', cdstr1]);

%ylabel(['days relative to Venus arrival date of ', cdstr2]);

%legend([h1(1), h2(1)], 'TOF (days)', 'DVT (km/sec)', 'Location','NorthWest');

% print -depsc -tiff -r300 dvt_tof.eps;