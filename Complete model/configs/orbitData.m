function orbit_data = orbitData()
    % Output a struct containing all the data relative to the spacecraft:
    %   - a:            Orbit's semi-major axis [km]
    %   - e:            Orbit's eccentricity [-]
    %   - i:            Orbit's inclination [rad]
    %   - T:            Orbit's period [s]
    %   - n_sc:         Average rotational rate of the orbit [rad/s]  
    %   - theta_G0:     Greenwich meridian initial latitude [rad]


    % Orbit semi-major axis [km]
    % a = 6716.488;
    a = 0.8016e4;
    % a = 16671;

    % Orbit eccentricity [-]
    e = 0.1678;
    % e = 0;
    % e = 0.00771;

    % Orbit inclination [rad]
    % i = pi/6;
    % i = deg2rad(50.3442);
    % i = deg2rad(40.901);
    i = 0;

    % Orbit's period [s]
    T = 2*pi*sqrt(a^3/astroConstants(13));

    % average rotational rate of the orbit
    n_sc = sqrt(astroConstants(13)/a^3);  
    % n_sc = 0.00114698;  

    % Greenwich meridian initial latitude [rad]
    theta_G0 = 0;

    orbit_data.a = a;
    orbit_data.e = e;
    orbit_data.i = i;
    orbit_data.T = T;
    orbit_data.n_sc = n_sc;
    orbit_data.theta_G0 = theta_G0;
end