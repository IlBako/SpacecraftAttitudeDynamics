function act_data = actuatorData()


    hr_max = 12.0; % [N*m*s]
    max_omega = 3600*2*pi/60; % [rad/s] 
    % hr_nom = hr_max/3; % [N*m*s] Nominal rotation speed
    hr_nom = 0;
    Ir = hr_max/max_omega;
    wr_nom = hr_nom/Ir; % [rad/s]

    hr_dot_max = 0.2; %[Nm]

    act_data.D_max = 80; %[Am^2] linear dipole moment 
    act_data.wr_nom = wr_nom; %[rad/s] nominal rotation speed
    act_data.hr_nom = hr_nom;
    act_data.hr_max = hr_max;
    act_data.hr_dot_max = hr_dot_max;
end