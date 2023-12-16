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

%% Run configs

configs;

%% Initial conditions after launcher release

% Initial angular velocity (body frame)
% in_cond.w0 = ((5-1)*rand(3,1) + 1) .* 1e-2;
in_cond.w0 = [3e-2 -1e-2 4e-2];

% Initial inertia wheel angular velocity
in_cond.wr0 = 0;

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
        actuator_data.max_dipole = 120; % [A*m^2]
        de_tumb = sim("Model.slx", sim_options);   

    case 2  % Pointing
        sim_options.StopTime = 'orbit_data.T';

    case 3  % De-tumbling + pointing
        algorithm = alg_vec{1};
        actuator_data.max_dipole = 120; % [A*m^2]
        de_tumb = sim("Model.slx", sim_options);

    case 4  % No Control
        algorithm = alg_vec{4};
        sim_options.StopTime = '10*orbit_data.T';
        actuator_data.max_dipole = 0; % [A*m^2]
        no_cont = sim("Model.slx", sim_options);

end

%% Plots

if strcmp(plot_gen, 'yes')
    generatePlots;
end