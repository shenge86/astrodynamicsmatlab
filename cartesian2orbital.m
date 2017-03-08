function [ORB] = cartesian2orbital(rovec, vovec)
% Converts Cartesian coordinates given by two arrays
% rovec - 3 dimensional position vector
% vovec - 3 dimensional velocity vector
% into the orbital elements
%
% Output vector ORB contains:
% a - semimajor axis (km)
% e - eccentricity (dimensionless)
% inc - inclination (radians)
% OMEGA - (radians)
% omega - (radians)
% f - true anomaly (radians)
% M - mean anomaly (radians) - optional


% Requires you to assign mu outside of the funtion
global mu

ro = norm(rovec);
vo = norm(vovec);

% a has to be negative
a = 1./(2./ro - vo.^2./mu);
if a < 0
    display('This orbit is not an ellipse.');
    return
end

hvec = cross(rovec,vovec);
cvec = cross(vovec, hvec) - mu.*rovec./ro;
e = norm(cvec)./mu;

% direction cosine matrix
% 3-1-3 Euler angles
% the 3 new axis
h = norm(hvec)
ie = cvec ./ (mu.*e) 
ih = hvec ./ h % normalized h vector
ip = cross(ih,ie)

norm(ie)
norm(ih)
norm(ip)

% part of a 3x3 matrix:
% C11 = ie(1);
% C12 = ie(2);
% C13 = ie(3);

% C21 = ip(1);
% C22 = ip(2);
% C23 = ip(3);

% C31 = ih(1);
% C32 = ih(2);
% C33 = ih(3);

OMEGA = atan2(ih(1),-ih(2))
inc = acos(ih(3))
omega = atan2(ie(3),ip(3))

sigma0 = ro.*vo./sqrt(mu)
% E = atan2(sigma0./sqrt(a),(1-ro./a))
E = acos((1-ro./a)./e)
f = 2.*atan(sqrt((1+e)/(1-e)).*tan(E./2))
% M = E - e.*sin(E) % - optional
ORB = [a e inc OMEGA omega f];