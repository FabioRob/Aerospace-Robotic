function [] = spacecraft_animation(assumption, planet_id_dep,planet_id_arr ,date_dep, date_arr, r0 ,v0 , animation_dir)

% -----------------------------------------------------------------------%
%	
% This function creates animation of spacecraft moving in transfer orbits
%
% Arguments:
%
% assumption      - 1/2/3 for real/circular&coplanar/coplanar orbits
% planet_id_dep   - ID of departure planet
% planet_id_arr   - ID of arrival planet
% date_dep        - departure date
% date_arr        - arrival date
% r0,v0           - state vector in transfer orbit at departure from planet_id_dep
% animation_dir   - directory in which the animation is stored
%
% -----------------------------------------------------------------------%

muSun=132712439935; % sun's gravitational parameter
    
JD_dep = date2JD(date_dep(1),date_dep(2),date_dep(3),date_dep(4),date_dep(5),date_dep(6));
JD_arr = date2JD(date_arr(1),date_arr(2),date_arr(3),date_arr(4),date_arr(5),date_arr(6));

tof=(JD_arr-JD_dep)*86400;

[rf,vf] = rv_from_r0v0( r0 ,v0,tof,muSun); % final spacecraft state vector
coe_fin = coe_from_sv(rf,vf,muSun);
TA_fin = coe_fin(6);

[J2000_coe, ~] = planetary_elements(assumption,max(planet_id_dep,planet_id_arr)) ;
a_outer = max(J2000_coe(1),coe_fin(7));

frame_rate=25;

if nargin==8
    mkdir(animation_dir);
    workingDir=animation_dir;
    h=figure('visible', 'off');
    k = 0;
    outputVideo = VideoWriter(fullfile(workingDir,[animation_dir,'movie.avi']));
    outputVideo.FrameRate = frame_rate;
    open(outputVideo)
end

hold all
whitebg('k')
view([0,90])

xlim([-1.5*a_outer 1.5*a_outer])
ylim([-1.5*a_outer 1.5*a_outer])

draw_sphere([0,0,0], 30, 'y', 1);

plotPlanetOrbit (3 , 2 , date_dep , 'g' , 1)
plotPlanetOrbit (3 , 3 , date_dep , 'w' , 1)
plotPlanetOrbit (3 , 4 , date_dep , 'r' , 1)
plotPlanetOrbit (3 , 5 , date_dep , 'b' , 1)
plotPlanetOrbit (3 , 6 , date_dep , 'y' , 1)
plotPlanetOrbit (3 , 7 , date_dep , 'c' , 1)

orbcol = [1,0.5,0] ;
OrbPlot(r0 ,v0,muSun,orbcol,1,TA_fin)

draw_sphere([0,0,0], 30, 'y', 1);

for JD = JD_dep:JD_arr
	
    t = (JD-JD_dep)*86400;
    [r,~] = rv_from_r0v0( r0 ,v0,t,muSun); % spacecraft position
    [year, month, day, hour, minute, second] = JD2date(JD);
    
    d = text(a_outer,a_outer,[num2str(day),'/',num2str(month),'/',num2str(year)]) ;
    
	pause(1/frame_rate);
    
    % plot spacecraft position
    spacecraft = draw_sphere(r, 10, 'w', 1);
    
    p=[];
    for i=2:7        
         [~, r_pl, ~, ~] = planet_elements_and_sv(3,i, year, month, day, hour, minute, second);
         p = [p,draw_sphere(r_pl, 10, 'g', 1)];       
    end
    
    if nargin==8      % create frames 
        filename=['frame',num2str(k,'%03d'),'.png'];
        frame = getframe(h);
        img = frame2im(frame);
        [imind,cm] = rgb2ind(img,256);
        fullname=[workingDir,'/',filename];
        imwrite(imind,cm,fullname);         
        
        writeVideo(outputVideo,img)
        
        if JD ~= JD_arr
            delete(spacecraft);
            delete(p);
            delete(d);
        end    
        
        k=k+1;
        
    end
    
end

% Create video(.avi) from frames
if nargin==8    
    close(outputVideo)
end

end



