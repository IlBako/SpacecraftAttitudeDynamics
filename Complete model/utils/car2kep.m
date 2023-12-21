function [a, e, i, OM, om, theta] = car2kep(rr, vv, mu)

% Transformation from Cartesian state to orbital elements
%
% [a, e, i, OM, om, theta] = rv2parorb (rr, vv, mu)
%
% Input arguments:
% ----------------------------------------------------------------
% rr            [3x1]   position vector                 [km]
% vv            [3x1]   velocity vector                 [km/s]
% mu            [1x1]   gravitational parameter         [km^3/s^2]
% 
% Output arguments:
% -----------------------------------------------------------------
% a             [1x1]   semi-major axis                 [km]
% e             [1x1]   eccentricity                    [-]
% i             [1x1]   inclination                     [rad]
% OM            [1x1]   RAAN                            [rad]
% om            [1x1]   pericenter anomaly              [rad]
% theta         [1x1]   true anomaly                    [rad]

I = [1 0 0];
J = [0 1 0];
K = [0 0 1];

r = norm(rr);                                       % Scalar distance 
v = norm(vv);                                       % Scalar Velocity

vr = dot(rr, vv)/r;                                 % Radial velocity

E = v^2/2 - mu/r;                                   % Specific mechanical energy
a = -mu/(2*E);                                      % Semi-major axis

h = cross(rr,vv);                                   % Specific angular momentum

i = wrapTo2Pi(acos(h(3)/norm(h)));                  % Inclination

N = cross(K, h);                                    % Nodes line

OM_check = i;
while OM_check > pi
    OM_check = OM_check - pi;
end

if OM_check== 0 
    OM = 0;
else
    if (N(2) < 0)
        OM = wrapTo2Pi((2*pi) - acos(N(1)/norm(N)));
    else
        OM = wrapTo2Pi(acos(N(1)/norm(N)));
    end
end

ee = cross(vv, h)/mu - rr/r;
e = norm(ee);

if OM_check == 0 || e < 1e-10       % APPROXIMATION HEEEEEEELP
    om = 0;
else
    if e < 1e-10
        om = 0;
    else
        if (ee(3) < 0)
            om = wrapTo2Pi(2*pi - acos(dot(N,ee)/(norm(N)*e)));
        else
            om = wrapTo2Pi(acos(dot(N,ee)/(norm(N)*e)));
        end
    end
end

if (vr < 0)
    theta = wrapTo2Pi(2*pi - real(acos(dot(ee,rr)/(r*e))));
else
    theta = wrapTo2Pi(acos(dot(ee,rr)/(r*e)));
end

end