%% SCRIPT CALCULATOR %%
% Calculates orbital elements eccentric anomaly and mean anomaly with respect to
% time

global mu
mu = 398600;

% initial conditions
t0 = 0;
E0 = 0;
M0 = 0;

rvec = [7100 0 1300];
vvec = [0 7.35 1];

ORB = cartesian2orbital(rvec,vvec)
a = ORB(1)
e = ORB(2)
inc = ORB(3)
OMEGA = ORB(4)
omega = ORB(5)
f = ORB(6)

% mean angular motion
n = sqrt(mu./a.^3);
period = 2.*pi./n; 
tf = period; % final time is the same as the period

% define time vector
sz = 1000; % number of elements
t = linspace(t0, tf, sz); 
M = zeros(1,sz);
E = zeros(1,sz);
E(1) = E0;
M(1) = M0;


for i=2:length(t)
    M(i) = n.*(t(i)-t(1));
    E(i) = kepler_E(e,M(i));
end

% convert to true anomaly
% f = atan(sqrt((1+e)./(1-e)).*tan(E./2)); % this will only cover from -90
% to 90 degrees
% f = 2.*atan(sqrt((1+e)./(1-e)).*tan(E./2));
f = 2.*atan2(sqrt(1-e).*cos(E./2), sqrt(1+e).*sin(E./2));
r = zeros(length(f),3);
v = zeros(length(f),3);
for i=1:length(t)
    [rovec, vovec] = orbital2cartesian(a, e, inc, OMEGA, omega, f(i));
    r(i,:) = rovec;
    v(i,:) = vovec;
end

% plot orbital trajectories in x, y, and z coordinates
grid on
plot3(r(:,1),r(:,2),r(:,3))
figure
plot(t,v)
legend('x-vel','y-vel','z-vel')
xlabel('Time (s)')
ylabel('Velociy (m/s)')

