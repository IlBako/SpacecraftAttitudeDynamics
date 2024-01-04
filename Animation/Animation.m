%% Animation

clc
close all

if ~exist('out','var')
    error("Run simulation first!")
end
warning("Satellite is not accurate for inertia and size, use for reference only")

A1 = out.A_BN;                  % first A matrix
A2 = out.A_BN_Sens;             % second A matrix
stlName = "Satellite.STL";      % stl file

percent_start = 40;         % simulation perc. for animation start. point
step_size = 100;
animation_length = 30000;    % number of frames displayed (total)

attitudeAnimation(A1, A2, stlName, percent_start, step_size, animation_length);


%% complete attitude animation

clc
close all

if ~exist('out','var')
    error("Run simulation first!")
end
warning("Satellite is not accurate for inertia and size, use for reference only")

A = out.A_BN_Sens;              % attitude
P = out.kepler_r_vec;           % radius of orbit
T = out.kepler_theta;           % theta of orbit
stlName = "Satellite.STL";

percent_start = 0;              % simulation perc. for animation start. point
step_size = 50;
animation_length = 0;        % number of frames displayed (total)

attitudeAnimation2(A, P, T, stlName, percent_start, step_size, animation_length, 1);


%% complete attitude animation - 3

clc
close all

if ~exist('out','var')
    error("Run simulation first!")
end
warning("Satellite is not accurate for inertia and size, use for reference only")

A = out.A_BN_Sens;              % attitude
P = out.kepler_r_vec;           % radius of orbit
T = out.kepler_theta;           % theta of orbit
MAG = out.sens_magn_B;
EHS = out.sens_ehs_r;
stlName = "Satellite.STL";

percent_start = 0;              % simulation perc. for animation start. point
step_size = 50;
animation_length = 5000;        % number of frames displayed (total)

attitudeAnimation3(A, P, T, MAG, EHS, stlName, percent_start, step_size, animation_length, 2);

