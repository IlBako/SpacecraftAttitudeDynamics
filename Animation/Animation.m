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

