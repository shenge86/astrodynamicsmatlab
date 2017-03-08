function [v, T] = orbCircular(r, M)
%% Calculates circular orbit's velocity and period given mass of planet and orbital radius
% Inputs:
% M - mass of planet (if not inputted then use Earth)
% r - orbital radius 
% Outputs:
% v - orbtial velocity (m/s)
% T - orbital period (minutes)

global G
G = 6.67e-11;

% default mass used if no M argument
if ~exist('M', 'var')
    M = 5.98e24; % mass of Earth (kg)
end

v = sqrt(G*M/r); % meters / second
T = 2*pi*r / v; % seconds
T = T/60; % minutes