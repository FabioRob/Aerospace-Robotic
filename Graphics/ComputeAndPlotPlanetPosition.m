function [coe, r, v, jd] = ComputeAndPlotPlanetPosition(assumption,planet_id,date,mar,col)

% -----------------------------------------------------------------------%
%
% ComputeAndPlotPlanetPosition computes the orbital elements,
% the state vector, and the julian day number at the specified date.
%
% Arguments :
%
% assumption - 1,2,3 for,repectively, real,coplanar and circular, coplanar orbits
% planet_id  - integer from 1 to 7 (from Mercury to Uranus)
% date       - [year,month,day,hour,minute,second]
% mar        - marker spec
% col        - color spec
%
% Outputs :
%
% coe        - [h e RA incl w TA a w_hat L M E] (with angles in degrees)
% r,v        - state vector
% jd         - julian day number
%
% -----------------------------------------------------------------------%

 year=date(1,1);
 month=date(1,2);
 day=date(1,3);
 hour=date(1,4);
 minute=date(1,5);
 second=date(1,6);
 format long
 [coe, r, v, jd]= planet_elements_and_sv(assumption,planet_id, year, month, day, hour, minute, second); % all angles in radians in coe    
    
 if nargin==5
       plot3(r(1),r(2),r(3),'Marker',mar,'color',col)
 elseif nargin==4
       plot3(r(1),r(2),r(3),'Marker',mar,'color','k')
 elseif nargin==3
       plot3(r(1),r(2),r(3),'Marker','o','color','k')
 end

end

