function [rovec, vovec] = orbital2cartesian(a, e, inc, OMEGA, omega, f)
% Converts the six orbital elements
%
% a - semimajor axis (km)
% e - eccentricity (dimensionless)
% inc - inclination (radians)
% OMEGA - (radians)
% omega - (radians)
% f - true anomaly (radians)
%
% into position and velocity vectors
% 
% rovec - 3 dimensional position vector (km)
% vovec - 3 dimensional velocity vector (km)

% Requires you to assign mu outside of the funtion

global mu

theta = omega + f;
p = a.*(1-e.^2);

r = p./(1+e.*cos(f)); % total distance (norm of r vector)
% v = sqrt(mu.*(2./r - 1./a));
h = sqrt(mu.*p); % magnitude of angular momentum

% inertial vector components of position and velocity are given by 
% 3-1-3 Euler angle set (OMEGA, inc, theta) applied.
rovec = r.*[cos(OMEGA).*cos(theta)-sin(OMEGA).*sin(theta).*cos(inc); sin(OMEGA).*cos(theta)+cos(OMEGA).*sin(theta).*cos(inc); sin(theta).*sin(inc)];
vovec = -mu./h .* [cos(OMEGA).*(sin(theta)+e.*sin(omega))+sin(OMEGA).*(cos(theta)+e.*cos(omega)).*cos(inc); sin(OMEGA).*(sin(theta)+e.*sin(omega))-cos(OMEGA).*(cos(theta)+e.*cos(omega)).*cos(inc); -(cos(theta)+e.*cos(omega)).*sin(inc)];

rovec = rovec';
vovec = vovec';