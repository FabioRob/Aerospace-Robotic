function [R] = R_to_perifocal(OM,i,om)

% -----------------------------------------------------------------------%
%
% R_to_perifocal function computes rotation matrix from planetocentri/heliocentric 
% frame to the perifocal frame.
%
% Arguments :
%
% OM  - right ascension of the ascending node(rad)
% i   - inclination of the orbit plane(rad)
% om  - argument of periapsis(rad)
%
% Output :
%
% R   - rotation matrix(3x3)
% 
% -----------------------------------------------------------------------%

R=[cos(OM)*cos(om)-sin(OM)*sin(om)*cos(i)  sin(OM)*cos(om)+cos(OM)*cos(i)*sin(om)  sin(i)*sin(om);
    -cos(OM)*sin(om)-sin(OM)*cos(i)*cos(om)  -sin(OM)*sin(om)+cos(OM)*cos(i)*cos(om)  sin(i)*cos(om);
    sin(OM)*sin(i)  -cos(OM)*sin(i)  cos(i)];

end

