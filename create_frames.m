function [kf] = create_frames (outputVideo , h ,frame_rate,JD_dep , JD_arr , r0 ,v0, workingDir , k  )

muSun=132712439935; % sun's gravitational parameter

for JD = JD_dep:JD_arr
	
    t = (JD-JD_dep)*86400;
    [r,~] = rv_from_r0v0( r0 ,v0,t,muSun); % spacecraft position
    [year, month, day, hour, minute, second] = JD2date(JD);    
        
	pause(1/frame_rate);
    
    % plot spacecraft position
    spacecraft = draw_sphere(r, 40, 'w', 1);
    
    p=[];
    for i=2:7        
         [~, r_pl, ~, ~] = planet_elements_and_sv(3,i, year, month, day, hour, minute, second);
         p = [p,draw_sphere(r_pl, 10, 'g', 1)];       
    end
    
      % create frames 
    filename=['frame',num2str(k,'%03d'),'.png'];
    frame = getframe(h);
    img = frame2im(frame);
    [imind,cm] = rgb2ind(img,256);
    fullname=[workingDir,'/',filename];
    imwrite(imind,cm,fullname);         
        
    writeVideo(outputVideo,img)
 
    delete(spacecraft);            
    delete(p);
        
    k=k+1;

end

kf=k ;

