%% Test function for cartesian to orbital and vice versa
clear
global mu
% deg = pi/180;
mu = 398600;
%%

% Input data
% r = [-6045 -3490 2500];
% v = [-3.457 6.618 2.533];

rvec = [7100 0 1300];
vvec = [0 7.35 1];

ORB = cartesian2orbital(rvec,vvec)
a = ORB(1)
e = ORB(2)
inc = ORB(3)
OMEGA = ORB(4)
omega = ORB(5)
f = ORB(6)

P = 2.*pi.*sqrt(a.^3./mu) % period of spacecraft
% [r2, v2] = orbital2cartesian(a, e, inc, OMEGA, omega, f)

% Propagate to future time
% numerically propagate set of 3 coupled differential equations
% r = norm(rvec)
% xdoubledot = - (mu ./ r^3) .* x
% similar for y and z

% t ranges from 0 to the period (loops around once)
% initial conditions for position are given by rvec
[t,x] = ode45(@orbfcn, [0, P], [rvec, vvec]);
xvec = x(:,1);
yvec = x(:,2);
zvec = x(:,3);
xdotvec = x(:,4);
ydotvec = x(:,5);
zdotvec = x(:,6);

plot3(x(:,1),x(:,2),x(:,3))
grid on
title('Orbit Curve')

