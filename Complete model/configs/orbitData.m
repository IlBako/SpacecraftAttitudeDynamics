function orbit_data = orbitData()
    % Output a struct containing all the data relative to the spacecraft:
    %   - a:     Orbit's semi-major axis [km]
    %   - e:     Orbit's eccentricity [-]
    %   - i:     Orbit's inclination [rad]
    %   - T:     Orbit's period [s]
    %   - n_sc:  Average rotational rate of the orbit [rad/s]  


    % Orbit semi-major axis [km]
    a = 3e4;

    % Orbit eccentricity [-]
    e = 0.3;

    % Orbit inclination [rad]
    i = pi/6;

    % Orbit's period [s]
    T = 2*pi*sqrt(a^3/astroConstants(13));

    % average rotational rate of the orbit
    n_sc = sqrt(astroConstants(13)/a^3);              

    orbit_data.a = a;
    orbit_data.e = e;
    orbit_data.i = i;
    orbit_data.T = T;
    orbit_data.n_sc = n_sc;
end