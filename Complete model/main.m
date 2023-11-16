clc; clearvars; close all;

sim_options.SolverType = 'Fixed-step';      % Set the solver type to Fixed-step
sim_options.Solver = 'ode4';                % Select ode4 as solver
sim_options.FixedStep = '0.1';              % Select a time step
sim_options.StartTime = '0';                % Start from 0 seconds [default]
% sim_options.StopTime = '5*9.952014160347722e+03';       % End time in seconds
sim_options.StopTime = '9.952014160347722e+03';       % End time in seconds

addpath("configs\");
addpath("utils\");

%% Data

astro_data = astronomicData;
sc_data = spaceCraftData;
orbit_data = orbitData;
in_cond = initialConditions(astro_data.n_eth);
pert_data = perturbationData(orbit_data.n_sc);

%% Simulation

out = sim("Model.slx", sim_options);

%% Post processing

figure();
plot(out.time, out.dynamics_omega, 'LineWidth', 1);
grid on;
xlabel("Time\ [s]", 'Interpreter','latex');
ylabel("$\omega$ [rad/s]", 'Interpreter','latex');
legend("$\omega_x$", "$\omega_y$", "$\omega_z$", 'Location','northeast', 'Interpreter', 'latex');

% tic
% A_det = zeros(1, length(out.time));
% for i = 1:length(out.time)
%     A_det(i) = det(out.kinematics_At(:, :, i));
% end
% disp(toc);
% figure
% plot(out.time, A_det);
% grid on;

figure
plot3(out.kepler_r_vec(:,1), out.kepler_r_vec(:,2), out.kepler_r_vec(:,3));
grid on; axis equal;