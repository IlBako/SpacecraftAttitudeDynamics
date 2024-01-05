clc; clearvars; close all;

addpath("Subsystems\");
addpath("configs\");
addpath("utils\");

%% Input handling

% Algorithm choice

alg_vec = {'De-tumbling', 'Pointing', 'De-tumbling + pointing', 'No Control'};

alg_in = 'Select the algorithm to use by typing the corresponding number: \n';
for i=1:length(alg_vec)
    alg_in = strcat(alg_in, num2str(i) + ". " + alg_vec{i} + " \n");
end
alg_idx = input(alg_in, "s");
alg_idx = str2double(alg_idx);

if alg_idx < 1 || alg_idx > length(alg_vec)
    error("Please insert a number between 1 and %d", length(alg_vec));
end

% Plots choice
plot_gen = input("Do you want to generate plots? (please answer with 'yes' or 'no'):  ", 's');

if strcmpi(plot_gen, 'yes')
    save_plots = input("Do you want to save the plots in either png or pdf?\n" + ...
        "(please answer with 'png', 'pdf' or 'no'):  ", 's');
end

% alg_idx = 4; plot_gen = 'yes'; save_plots = 'no';

%% Run configs

configs;

%% Initial conditions after launcher release

% Initial angular velocity (body frame)
% in_cond.w0 = ((5-1)*rand(3,1) + 1) .* 1e-2;
% in_cond.w0 = [3e-2 -1e-2 4e-2];
in_cond.w0 = [pi/3 pi/6 pi/4.5];

% Initial inertia wheel angular velocity
in_cond.wr0 = 0;

in_cond.q0 = dcm2quat(in_cond.A0);

pointing_toll = deg2rad(5);
pointing_k1 = 1*max(diag(sc_data.I_mat))/min(diag(sc_data.I_mat))*ones(1,3);
pointing_k2 = -1e-3*diag(sc_data.I_mat) .* [1 1 1]';
% pointing_k2 = 0;
% pointing_k1 = 1*max(diag(sc_data.I_mat))/min(diag(sc_data.I_mat))*ones(1,3);
% pointing_k2 = -1e-3*diag(sc_data.I_mat) .* [1 1 1]';
% pointing_k1 = sc_data.I_mat(3,3)/sc_data.I_mat(1,1)*[1 1 1];
% pointing_k2 = -0.5*diag(sc_data.I_mat) .* [1 1 1]';

[rr, vv] = kep2car(orbit_data.a, orbit_data.e, orbit_data.i, 0,0,0, astro_data.muE);
in_z = rr/norm(rr); in_y = cross(rr, vv)/norm(cross(rr, vv)); in_x = cross(in_y, in_z)/norm(cross(in_y, in_z));
A_LN0 = [in_x' in_y' in_z'];

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
        de_tumb = sim("Model.slx", sim_options);   

    case 2  % Pointing
        algorithm = alg_vec{2};
        sim_options.StopTime = '25*orbit_data.T';
        in_cond.w0 = [-2e-3 3e-3 -5e-4];
        % in_cond.A0 = A_LN0;
        point = sim("Model.slx", sim_options);

    case 3  % De-tumbling + pointing
        algorithm = alg_vec{1};
        actuator_data.hr_nom = 0;
        de_tumb = sim("Model.slx", sim_options);

    case 4  % No Control
        algorithm = alg_vec{4};
        sim_options.StopTime = '5*orbit_data.T';
        in_cond.w0 = [4e-4 3e-3 3e-4];
        no_cont = sim("Model.slx", sim_options);

end

%% Plots

if strcmp(plot_gen, 'yes')
    generatePlots;
end