function E = kepler_E(e,M)
% Uses Newton's method to solve Kepler's equation for an elliptical orbit
% given eccentricity and mean anomaly.
% Mean anomaly will need to be derived from sqrt(mu/a^3)*(t-t0)
% which is outside the scope of this function.
%
%

% error tolerance
error = 1.e-8;

% select starting value for E
if M < pi
    E = M + e/2;
else
    E = M - e/2;
end

% iterate until E is determined within error tolerance
% note 1 - e*cos(E) is the derivate of M = E-e*sinE
ratio = 1;
while abs(ratio) > error
    ratio = (E-e*sin(E) - M) / (1 - e*cos(E));
    E = E - ratio;
end