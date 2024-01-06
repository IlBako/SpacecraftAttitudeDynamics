clc; clearvars; close all;

addpath("Subsystems\");
addpath("configs\");
addpath("utils\");

%% Input handling

% Algorithm choice
alg_vec = {'De-tumbling', 'Pointing', 'Detumbling and Pointing', 'No Control'};

alg_in = 'Select the algorithm to use by typing the corresponding number: \n';
for i=1:length(alg_vec)
    alg_in = strcat(alg_in, num2str(i) + ". " + alg_vec{i} + " \n");
end
alg_idx = input(alg_in, "s");
alg_idx = str2double(alg_idx);

if alg_idx < 1 || alg_idx > length(alg_vec) || isnan(alg_idx)
    error("Please insert a number between 1 and %d", length(alg_vec));
end

% Plots choice
keepAsking1 = 1;
save_plots = "no";
while keepAsking1
    plot_gen = input("Do you want to generate plots? (please answer with 'yes' or 'no'):  ", 's');    
    if strcmpi(plot_gen, 'yes') || strcmpi(plot_gen, 'y')
        keepAsking1 = 0;
    elseif strcmpi(plot_gen, 'no') || strcmpi(plot_gen, 'n')
        keepAsking1 = 0;
    end
end

clc
fprintf("Algorithm: '" + alg_vec(alg_idx) + "', plots generation: '" + plot_gen + "'\n")
disp("Running simulation...")
tic

%% Run configs

configs;

%% Initial conditions after launcher release

% Initial angular velocity (body frame)
in_cond.w0 = deg2rad([7.5 7.5 7.5]);

% Initial dcm attitude [-]
in_cond.A0 = eye(3);
% Initial position on orbit [rad]
in_cond.theta0 = 0;

% Initial inertia wheel angular velocity
in_cond.wr0 = 0;

[in_cond.r0, in_cond.v0] = kep2car(orbit_data.a, orbit_data.e, orbit_data.i, 0,0,0, astroConstants(13));

%% Control logic parameters
pointing_toll = deg2rad(5);
pointing_k1 = 1*max(diag(sc_data.I_mat))/min(diag(sc_data.I_mat))*ones(1,3);
pointing_k2 = -1e-3*diag(sc_data.I_mat) .* [1 1 1]';
point_stop_time = 2*orbit_data.T;

%% Setup simulink options

% Simulink options
sim_options.SolverType = 'Fixed-step';      % Set the solver type to Fixed-step
sim_options.Solver = 'ode4';                % Select ode4 as solver
sim_options.FixedStep = '0.1';              % Select a time step less than or equal to the minimum step size
sim_options.StartTime = '0';                % Start from 0 seconds [default]
sim_options.StopTime = '100*orbit_data.T';      % End time in seconds

%% Run simulation

switch alg_idx

    case 1  % De-tumbling
        algorithm = alg_vec{1};
        actuator_data.hr_nom = 0;
        detumb = sim("Model.slx", sim_options);   

    case 2  % Pointing
        algorithm = alg_vec{2};
        sim_options.StopTime = '15*orbit_data.T';
        in_cond.w0 = [5.5197e-04 0.0063 0.0078];
        point = sim("Model.slx", sim_options);

    case 3 % Detumbling and pointing
        algorithm = alg_vec{1};
        actuator_data.hr_nom = 0;
        detumb = sim("Model.slx", sim_options);
        disp("Detumbling finished in " + toc + "seconds");

        algorithm = alg_vec{2};
        in_cond.w0 = detumb.dynamics_omega(end,:);
        orbit_data.theta_G0 = detumb.theta_G(end);
        in_cond.A0 = detumb.A_BN(:,:,end);
        in_cond.theta0 = detumb.kepler_theta(end);
        in_cond.r0 = detumb.kepler_r_int(end,:);
        in_cond.v0 = detumb.kepler_vel(end,:);
        point_stop_time = 2*orbit_data.T;

        sim_options.StopTime = '15*orbit_data.T';
        point = sim("Model.slx", sim_options);

    case 4  % No Control
        algorithm = alg_vec{4};
        sim_options.StopTime = '1*orbit_data.T';
        in_cond.w0 = [2.2951e-04 8.7378e-03 -2.6631e-04];
        no_cont = sim("Model.slx", sim_options);

end

%% Plots

if strcmp(plot_gen, 'yes') || strcmp(plot_gen, 'y')
    generatePlots;
end

disp("Simulation finished! Time elapsed: " + toc + " seconds")

%% Animation

switch alg_idx
    case 1  % De-tumbling
        out = detumb;   
    case 2  % Pointing
        out = point;
    case 3  % De-tumbling + pointing
        out = point;
    case 4  % No Control
        out = no_cont;
end

percent_start = 90;
step_size = 150;
animation_length = 0;
cam_choice = 1;

attitudeAnimation(out.A_BN_Sens, out.kepler_r_int, out.kepler_theta, "Satellite.STL", percent_start, step_size, animation_length, cam_choice);

