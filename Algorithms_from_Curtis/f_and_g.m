function [f, g] = f_and_g(x, t, ro, a,mu)
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
%
% This function calculates the Lagrange f and g coefficients.
%
% mu - the gravitational parameter (kmˆ3/sˆ2)
% a - reciprocal of the semimajor axis (1/km)
% ro - the radial position at time t (km)
% t - the time elapsed since t (s)
% x - the universal anomaly after time t (kmˆ0.5)
% f - the Lagrange f coefficient (dimensionless)
% g - the Lagrange g coefficient (s)
%
% User M-functions required: stumpC, stumpS
% ------------------------------------------------------------

z = a*x^2;
%...Equation 3.66a:
f = 1 - x^2/ro*stumpC(z);
%...Equation 3.66b:
g = t - 1/sqrt(mu)*x^3*stumpS(z);