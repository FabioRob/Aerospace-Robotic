function []  = plotApseLine(oe,mu,col)
% This function plots the apse line of the orbit 
%  oe - orbital elements:
%  [h e RA incl w TA]
% where
% h = angular momentum (kmˆ2/s)
% e = eccentricity
% RA = right ascension of the ascending node (rad)
% incl = inclination of the orbit (rad)
% w = argument of perigee (rad)
% TA = true anomaly (rad)
% mu - gravitational parameter of the central body 

h = oe(1);
e=oe(2);
RA=oe(3);
incl=oe(4);
w=oe(5);
TA=oe(6); % not used

% compute periapsis and apoapsis radii in perifocal frame
rp_per = (h^2 /  (mu*(1+e))) * [1,0,0]';
ra_per = (h^2 /  (mu*(1-e))) * [-1,0,0]';

% compute periapsis and apoapsis radii in central body frame
Rot_to_centralFrame = (R_to_perifocal(RA,incl,w))' ;
rp_central = Rot_to_centralFrame * rp_per;
ra_central = Rot_to_centralFrame * ra_per;
apseLine=[rp_central';ra_central'];
hold on
plot3(apseLine(:,1),apseLine(:,2),apseLine(:,3),'--','color',col)
plot3(rp_central(1),rp_central(2),rp_central(3),'+','color',col) % perihelium


