function act_data = actuatorData()

    hr_max_point = 12; % [N*m*s]
    hr_max_slew = 12/6; % [N*m*s]
    max_omega = 4000*2*pi/60; % [rad/s] 
    hr_nom = hr_max_point/3; % [N*m*s] Nominal rotation speed
    Ir = hr_max_point/max_omega;
    wr_nom = hr_nom/Ir; % [rad/s]

    hr_dot_max = 0.2; %[Nm]

    act_data.D_max = 140; %[Am^2] linear dipole moment 
    act_data.wr_nom = wr_nom; %[rad/s] nominal rotation speed
    act_data.hr_nom = hr_nom;
    act_data.hr_max_point = hr_max_point;
    act_data.hr_max_slew = hr_max_slew;
    act_data.hr_dot_max = hr_dot_max;
end