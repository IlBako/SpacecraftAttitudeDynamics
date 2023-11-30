function astro_data = astronomicData()
    % Output a struct containing all the astronomic data useful for the simulations:
    %   - muE:      Planetary constants of Earth (mu = G*mass) [km^3/s^2]
    %   - c:        Speed of light in the vacuum [km/s]
    %   - R_earth:  Earth's mean radius [Km]
    %   - n_eth:    Average rotational rate using Earth's radius [rad/s]
    %   - n_sun:    Average rotational rate of Earth's orbit around the Sun [rad/s]
    %   - eps:      Earth's orbit inclination [rad]
    %   - R_sun:    Earth's orbit mean radius (1AU) [km]
    %   - wE:       Earth's angular velocity around it's axis [rad/s]

    % Planetary constants of Earth (mu = G*mass) [km^3/s^2]
    muE = astroConstants(13);

    % Speed of light in the vacuum [m/s]
    c = astroConstants(5)*1e3;

    % Earth's mean radius [Km]
    R_earth = astroConstants(23);

    % Average rotational rate using Earth's radius [rad/s]
    % n_eth = 15.04*pi/180*1/3600; % DA RIVEDERE STO MAGIC NUMBER POPI POPI
    % n_eth = 2*pi/(24*3600);
    n_eth = sqrt(muE/R_earth^3);

    % Average rotational rate of Earth's orbit around the Sun [rad/s]
    n_sun = 2*pi/(365*24*60^2);

    % Earth's orbit inclination [rad]
    epsilon = deg2rad(23.45);

    % Earth's orbit mean radius (1AU) [km]
    r_sun = astroConstants(2);

    % Earth angular velocity around it's axis [rad/s]
    wE = deg2rad(15.04)/3600;

    astro_data.muE = muE;
    astro_data.c = c;
    astro_data.R_earth = R_earth;
    astro_data.n_eth = n_eth;
    astro_data.n_sun = n_sun;
    astro_data.eps = epsilon;
    astro_data.R_sun = r_sun;
    astro_data.wE = wE;
end