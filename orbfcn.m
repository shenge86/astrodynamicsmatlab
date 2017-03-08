function dx = orbfcn(t,x)
% OUTPUT:
% dx is left side of set of differential equations
% corresponding to derivatives  
%
% INPUT:
% equation on the right side is what each derivative is equal to

% Set of 3 coupled differential equations of second order
% decomposed into set of 6 coupled differential equations of first order
% inputs are as follows:
% x(1) - x
% x(2) - y
% x(3) - z
% x(4) - xdot
% x(5) - ydot
% x(6) - zdot
% outputs:
% dx(1) - xdot
% dx(2) - ydot
% dx(3) - zdot
% dx(4) - xdoubledot
% dx(5) - ydoubledot
% dx(6) - zdoubledot
global mu

dx(1) = x(4); % xdoubledot = d/dt (xdot)
dx(4) = -mu.*x(1)./norm([x(1) x(2) x(3)]);% xdoubledot = - (mu ./ r^3) .* x

dx(2) = x(5);
dx(5) = -mu.*x(2)./norm([x(1) x(2) x(3)]);

dx(3) = x(6);
dx(6) = -mu.*x(3)./norm([x(1) x(2) x(3)]);

dx = dx(:);