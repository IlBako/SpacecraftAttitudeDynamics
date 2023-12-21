clc, clearvars
close all

addpath("configs\");
addpath("utils\");
addpath("Subsystems\");

%% Load data

astro_data = astronomicData;
sc_data = spaceCraftData;
orbit_data = orbitData;
in_cond = initialConditions(astro_data.n_eth);
pert_data = perturbationData(orbit_data.n_sc);
sensor_data = sensorData;
worldMag_data = load("data\WMD.mat");

% For pointing n should be replaced by theta dot when calculating the LVLH
% frame

A_target=[1,0,0;0,1,0;0,0,1]; %vogliamo puntare la terra
algorithm = 'No control';
in_cond.q0 = dcm2quat(in_cond.A0);

% Kpx=-1;
% Kpy=-1;
% Kpz=-1;
% 
% Kdx=-1e1;
% Kdy=-1e1;
% Kdz=-1e1;

pointing_k1 = -sc_data.I_mat(3,3)/sc_data.I_mat(1,1)*[1 1 1];
pointing_k2 = -0.5*diag(sc_data.I_mat) .* [1 1 1]';

%% Simulation

% Simulink options
sim_options.SolverType = 'Fixed-step';      % Set the solver type to Fixed-step
sim_options.Solver = 'ode4';                % Select ode4 as solver
sim_options.FixedStep = '0.1';              % Select a time step less than or equal to the minimum step size
sim_options.StartTime = '0';                % Start from 0 seconds [default]
sim_options.StopTime = '1*orbit_data.T';      % End time in seconds

% Call simulink model
out = sim("Model.slx", sim_options);

%% Post processing

% Plot - om
figure()
plot(out.time, out.dynamics_omega, 'LineWidth', 1);
grid on;
title("Angular velocities");
xlabel("Time\ [s]", 'Interpreter','latex');
ylabel("$\omega$ [rad/s]", 'Interpreter','latex');
legend("$\omega_x$", "$\omega_y$", "$\omega_z$", 'Location','northeast', 'Interpreter', 'latex');

% Plot - disturbs
% figure()
% plot(out.time, vecnorm(out.perturbations_I_M_GG, 2, 2), 'LineWidth', 1);
% grid on, hold on
% plot(out.time, vecnorm(out.perturbations_I_M_Mag, 2, 2), 'LineWidth', 1);
% plot(out.time, vecnorm(out.perturbations_I_M_SRP, 2, 2), 'LineWidth', 1);
% xlabel("Time\ [s]", 'Interpreter','latex');
% ylabel("$\dot{omega}$ [rad/s]", 'Interpreter','latex');
% legend("$GG$", "$MAG$", "$SRP$", 'Location','northeast', 'Interpreter', 'latex');
% title("Disturbances");

% Plot - kinetic energy
figure()
T = 0.5 * (sc_data.I_mat(1,1)*out.dynamics_omega(:,1).^2 + sc_data.I_mat(2,2)*out.dynamics_omega(:,2).^2 + sc_data.I_mat(3,3)*out.dynamics_omega(:,3).^2);
plot(out.time, T, 'LineWidth', 1)
grid on
title("Kinetic energy")
xlabel("Time\ [s]", 'Interpreter','latex');
ylabel("T [J]", 'Interpreter','latex');
legend("Kinetic Energy", 'Location','northeast', 'Interpreter', 'latex');

% Plot - err om_BL
% figure
% plot(out.time, out.w_BL - out.dynamics_omega, 'LineWidth', 1) % TODO check if the formula is correct
% grid on
% title("Angular velocity error")
% xlabel("Time\ [s]", 'Interpreter','latex');
% ylabel("$\omega$ [rad/s]", 'Interpreter','latex');
% legend("$\omega_{x} error$", "$\omega_{y} error$", "$\omega_{z} error$", 'Location','northeast', 'Interpreter', 'latex');

% Plot - stability graph
% I_x = sc_data.I_mat(1,1);
% I_y = sc_data.I_mat(2,2);
% I_z = sc_data.I_mat(3,3);
% 
% K_x = (I_z - I_y)/I_x;
% K_y = (I_z - I_x)/I_y;
% K_z = (I_y - I_x)/I_z;
% 
% figure()
% plot(K_y, K_x, 'or')
% hold on
% plot([-1 1], [-1 1])
% plot([-1 1], [0 0], color = '[0.2, 0.2, 0.2]')
% plot([0 0], [-1 1], color = '[0.2, 0.2, 0.2]')
% grid on, axis equal
% xlim([-1 1])
% ylim([-1 1])
% title("Stability plot")

% Plot - orbit
figure()
earthPlot();
plot3(out.kepler_r_vec(:,1), out.kepler_r_vec(:,2), out.kepler_r_vec(:,3), 'LineWidth', 1);
grid on, axis equal
view(30, 20)
title("Orbit")
