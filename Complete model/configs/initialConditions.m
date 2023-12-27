function in_cond = initialConditions(n_eth)
    % Output a struct containing all the initial conditions:
    %   - w0:       Initial angular velocity [rad/s]
    %   - A0:       Initial dcm attitude [-]
    %   - theta0:   Initial position on orbit [rad]
    %   - wr0:      Initial inertia wheel angular velocity [rad/s]

    % Initial angular velocity [rad/s]
    % w0 = [1e-6 1e-6 n_eth];
    % w0 = [3e-2 1e-2 3e-2];
    % w0 = [pi/2 pi/3 pi/6];
    % w0 = [0.45 0.52 0.55];
    w0 = [0.5e-3 1e-4 0.9e-3];
    % Initial dcm attitude [-]
    A0 = eye(3);
    % Initial position on orbit [rad]
    theta0 = 0;
    % Initial inertia wheel angular velocity [rad/s]
    % wr0 = 0;
    wr0 = 8.8234;

    in_cond.w0 = w0;
    in_cond.A0 = A0;
    in_cond.theta0 = theta0;
    in_cond.wr0 = wr0;
end