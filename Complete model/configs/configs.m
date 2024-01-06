rng("default")

%% Astronomic data

% Planetary constants of Earth (mu = G*mass) [km^3/s^2]
astro_data.muE = astroConstants(13);
% Speed of light in the vacuum [m/s]
astro_data.c = astroConstants(5)*1e3;
% Earth's mean radius [Km]
astro_data.R_earth = astroConstants(23);
% Average rotational rate using Earth's radius [rad/s]
astro_data.n_eth = sqrt(astro_data.muE/astro_data.R_earth^3);
% Average rotational rate of Earth's orbit around the Sun [rad/s]
astro_data.n_sun = 2*pi/(365*24*60^2);
% Earth's orbit inclination [rad]
astro_data.eps = deg2rad(23.45);
% Earth's orbit mean radius (1AU) [km]
astro_data.R_sun = astroConstants(2);
% Earth angular velocity around it's axis [rad/s]
astro_data.wE = deg2rad(15.04)/3600;

%% Spacecraft data

% Principal moments of Intertia [Kg*m^2]
I = [23.4677 32.3556 32.9733];
% Principal inertia matrix [Kg*m^2]
sc_data.I_mat = I.*eye(3);
% Inverse of intertia matrix
sc_data.I_inv = inv(sc_data.I_mat);
% NOTE: The inverse is calculated before hand to save computation time
% as the inertia matrix does not change during the mission

% Normals to body 
sc_data.NB = [ 1  0  0;      % Body 1       
               0  1  0;      % Body 2
              -1  0  0;      % Body 3
               0 -1  0;      % Body 4
               0  0  1;      % Body 5
               0  0 -1;      % Body 6
               0  0  1;      % Panel 1
               0  0 -1;      % Panel 2
               0  0  1;      % Panel 3
               0  0 -1];     % Panel 4;

% specular reflection coefficient
rhoS_body = 0.5;
rhoS_panel = 0.1;
sc_data.rhoS = [repmat(rhoS_body, 6, 1); repmat(rhoS_panel, 4, 1)];
% diffusion reflection coefficient
sc_data.rhoD = repmat(0.1, 10, 1);
% Surface area [m^2]
A_body = [0.7088;0.7088;0.7088;0.7088;0.64;0.64];
A_panel = 96e-2*ones(4,1);
sc_data.A = [A_body; A_panel];

% distance from baricentre [m]
sc_data.rF = [ 40    0       0;      % Body 1
                0   40       0;      % Body 2
              -40    0       0;      % Body 3
                0  -40       0;      % Body 4
                0    0    44.3;      % Body 5
                0    0   -44.3;      % Body 6
                0    0     100;      % Panel 1
                0    0     100;      % Panel 2
                0    0    -100;      % Panel 3
                0    0    -100];     % Panel 4;
sc_data.rF = sc_data.rF * 1e-2;
% Number of faces
sc_data.NumFaces = size(sc_data.NB, 1);
% Drag coefficient
sc_data.cD = 2.1;

%% Orbit data

% Orbit semi-major axis [km]
orbit_data.a = 0.8016e4;
% Orbit eccentricity [-]
orbit_data.e = 0.1678;
% Orbit inclination [rad]
orbit_data.i = 0;
% Orbit's period [s]
orbit_data.T = 2*pi*sqrt(orbit_data.a^3/astroConstants(13));
% average rotational rate of the orbit
orbit_data.n_sc = sqrt(astroConstants(13)/orbit_data.a^3);
% Greenwich meridian initial latitude [rad]
orbit_data.theta_G0 = 0;

%% Perturbation data

% Solar radiation intensity [W/m^2]
pert_data.Fe = 1358;

% Density data -> rho = rho_0 * exp(-(h-h0)/H)
pert_data.h0_vect = [0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 ...
        180 200 250 300 350 400 450 500 600 700 800 900 1000]';
pert_data.rho0_vect = [1.225 3.899*1e-2 1.774*1e-2 3.972*1e-3 1.057*1e-3 ...
        3.206*1e-4 8.770*1e-5 1.905*1e-5 3.396*1e-6 5.297*1e-7 9.661*1e-8 ...
        2.438*1e-8 8.484*1e-9 3.845*1e-9 2.070*1e-9 5.464*1e-10 ...
        2.789*1e-10 7.248*1e-11 2.418*1e-11 9.158*1e-12 3.725*1e-12 ...
        1.585*1e-12 6.967*1e-13 1.454*1e-13 3.614*1e-14 1.170*1e-14 ...
        5.245*1e-15 3.019*1e-15]';
pert_data.H_vect = [7.249 6.349 6.682 7.554 8.382 7.714 6.549 5.799 5.382 5.877 ...
        7.263 9.473 12.636 16.149 22.523 29.740 37.105 45.546 53.628 ...
        53.298 58.515 60.828 63.822 71.835 88.667 124.64 181.05 268.00]';

% Residual dipole of the spacecraft [A*m^2]
pert_data.dip_sc = [0.05 0.05 0.05];

% Dipole model
g01 = -29404.8;
g11 = -1450.9;
h11 = 4652.5;
pert_data.H0 = sqrt(g01^2 + g11^2 + h11^2);

%% Sensor data

%%% Magnetometer

% [T] Full scale measurement of the sensor (range is +-)
magnetometer.FSS = 100e-6;
% [%] Accuracy of the sensor
magnetometer.acc = 0.05;
% [%] Non-linearity of the sensor
magnetometer.lin = 0.015;
% [V/T] Sensitivity of the sensor
magnetometer.sens = 100*1e3;
% [T] Standard deviation of the sensor
magnetometer.std_dev = sqrt(magnetometer.acc ^ 2 + magnetometer.lin^2)/sqrt(3) * magnetometer.FSS / 100;
% [T^2] Variance of the sensor
magnetometer.variance = magnetometer.std_dev^2;
% [V] Quantization interval in Volt
V_quant = 0.025*2;
% [T] Quantization in Tesla
magnetometer.T_quant = V_quant/magnetometer.sens;
% [Hz] Frequency of the sensor
magnetometer.freq = 10;

% To calculate the magnetometer non orthogonality error 
% we use the inverse approach of hard and soft iron calibration.
% The way hard and soft iron calibration work is by removing the offset
% of the center of the magnetometer sphere and by transforming the
% noisy magnetometer measure, which creates an ellipse, in a sphere

a = deg2rad(rand(3,1)*1); % Soft iron - orthogonality within 1 deg
b = deg2rad(rand(3,1)*360); % Hard iron - rotation over 360 deg
% Orthogonality error matrix: this matrix both rotates and stretches
% the magnetic field to turn it into an ellipse with offset center
magnetometer.A_ortho = [cos(a(1))               sin(a(1))*cos(b(1))     sin(a(1))*sin(b(1));
                        sin(a(2))*cos(b(2))     cos(a(2))               sin(a(2))*sin(b(2));
                        sin(a(3))*cos(b(3))     sin(a(3))*sin(b(3))     cos(a(3))];

%%% Earth horizon sensor
% [rad] Accuracy of the sensor (0.05deg at 3 sigma)
horizon.acc = deg2rad(0.05)/3;
% Misalignment error
horizon.misalign = deg2rad(0.5);
% [Hz] Frequency of the sensor
horizon.freq = 10;

sensor_data.magnetometer = magnetometer;
sensor_data.horizon = horizon;

%% Actuator data

%%% Magnetorquers
%[Am^2] linear dipole moment
actuator_data.D_max = 400;  

%%% Inertia wheel
% Maximum rotational velocity
max_omega = 3600*2*pi/60; % [rad/s] 
% Maximum hr for pointing
actuator_data.hr_max_point = 12;
% [N*m*s] Nominal rotation speed for no control
actuator_data.hr_nom = actuator_data.hr_max_point/3;
% Inertia of the wheel
Ir = actuator_data.hr_max_point/max_omega;
% Maximum hr for slew
actuator_data.hr_max_slew = 12/6;
% Maximum hr_dot
actuator_data.hr_dot_max = 0.2;
% [rad/s] nominal rotation speed
actuator_data.wr_nom = actuator_data.hr_nom/Ir;

%% Magnetic field 13th order data
worldMag_data = load("data\WMD.mat");