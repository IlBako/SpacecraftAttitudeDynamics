function orbit_data = orbitData()
    % Output a struct containing all the data relative to the spacecraft:
    %   - a:     Orbit's semi-major axis [km]
    %   - e:     Orbit's eccentricity [-]
    %   - i:     Orbit's inclination [rad]
    %   - n_sc:  Average rotational rate of the orbit [rad/s]  


    % Orbit semi-major axis [km]
    a = 1e4;

    % Orbit eccentricity [-]
    e = 0.5;

    % Orbit inclination [rad]
    i = pi/6;

    % average rotational rate of the orbit
    n_sc = sqrt(astroConstants(13)/a^3);              

    orbit_data.a = a;
    orbit_data.e = e;
    orbit_data.i = i;
    orbit_data.T = 2*pi*sqrt(a^3/astroConstants(13));
    orbit_data.n_sc = n_sc;
end