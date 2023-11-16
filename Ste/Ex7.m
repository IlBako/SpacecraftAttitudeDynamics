%% task 1 to 4

clear, clc
close all
% close_system

om_earth = 15.04 * (pi/180) * 1/3600;
R_earth = astroConstants(23);

e = 0.5;
a = 10000;
i = pi/6;
mu = astroConstants(13);
n = sqrt(mu/(astroConstants(23)^3));
n_orb = sqrt(mu/(a^3));
c = astroConstants(5);

Rs = astroConstants(3); % Sun radius

% phisical par
Nb_vers = [1 0 0; 0 1 0; -1 0 0; 0 -1 0; 0 0 1; 0 0 -1; 1 0 0; -1 0 0; 1 0 0; -1 0 0];
rho_s = [0.5 0.5 0.5 0.5 0.5 0.5 0.1 0.1 0.1 0.1]';
rho_d = 0.1 * ones(10,1);
A = 0.01 .* [6 6 6 6 4 4 12 12 12 12]';
rf = [10 0 0; 0 10 0; -10 0 0; 0 -10 0; 0 0 15; 0 0 -15; 0 0 45; 0 0 45; 0 0 -45; 0 0 -45] .* 1e-2;
jB = [0.01 0.05 0.01]' .* 1e-9;

I_x = 100.9 * 1e-2;
I_y = 25.1 * 1e-2;
I_z = 91.6 * 1e-2;

om_x0 = 1e-6;
om_y0 = 1e-6;
om_z0 = n;

om_x_L = 0;
om_y_L = 0;
om_z_L = n_orb;
om_L = [om_x_L om_y_L om_z_L]';

A0 = eye(3);
om0 = [om_x0 om_y0 om_z0];
J = [I_x 0 0; 0 I_y 0; 0 0 I_z];
J_inv = inv(J);

T_orbit = 2*pi*sqrt((a^3)/mu);

sim_options.Solver = 'ode4';                % Select ode4 as solver
sim_options.SolverType = 'Fixed-step';      % Set the solver type to Fixed-step
sim_options.FixedStep = '0.5';              % Select a time step of 0.1 s
sim_options.StartTime = '0';                % Start from 0 seconds [default]
sim_options.StopTime = 'T_orbit';           % End the simulation at t=10s

out=sim('SimEx1.slx', sim_options);


%% figure

plot(out.time, out.om) 
grid on
legend("om x", "om y", "om z", Location="best")

T = 0.5 * (I_x*out.om(:,1).^2 + I_y*out.om(:,2).^2 + I_z*out.om(:,3).^2);
figure
plot(out.time, T)
grid on
legend("T", Location="best")

figure
plot(out.time, out.err_om_BL)
grid on
legend("err om x", "err om y", "err om z", Location="best")

K_x = (I_z - I_y)/I_x;
K_y = (I_z - I_x)/I_y;
K_z = (I_y - I_x)/I_z;

figure
plot(K_y, K_x, 'or')
hold on
plot([-1 1], [-1 1])
plot([-1 1], [0 0], color = '[0.2, 0.2, 0.2]')
plot([0 0], [-1 1], color = '[0.2, 0.2, 0.2]')
grid on, axis equal
xlim([-1 1])
ylim([-1 1])


%% orbit plot

plot3(out.R(:,1), out.R(:,2), out.R(:,3))
axis equal, grid on

