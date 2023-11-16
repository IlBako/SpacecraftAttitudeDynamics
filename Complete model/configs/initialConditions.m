function in_cond = initialConditions(n_eth)
    % Output a struct containing all the initial conditions:
    %   - w0:       Initial angular velocity [rad/s]
    %   - A0:       Initial dcm attitude [-]
    %   - theta0:   Initial position on orbit [rad]

    % Initial angular velocity [rad/s]
    w0 = [1e-6 1e-6 n_eth];
    % w0 = [0.45 0.52 0.55];
    % Initial dcm attitude [-]
    A0 = eye(3);
    % Initial position on orbit [rad]
    theta0 = 0;

    in_cond.w0 = w0;
    in_cond.A0 = A0;
    in_cond.theta0 = theta0;
end