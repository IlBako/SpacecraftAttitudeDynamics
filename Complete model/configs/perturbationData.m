function pert_data = perturbationData()  
    % Output a struct containing all the data relative to the perturbations:
    % GG Perturbation:
    %   - w_LN:     Reference moving angular velocity [rad/s]
    % SRP Perturbation
    %   - Fe:       Solar radiation intensity [W/m^2]
    % Magnetic Perturbation:
    %   - H0:        
    %   - m:

    % Solar radiation intensity [W/m^2]
    Fe = 1358;

    %% Air drag

    % Density data -> rho = rho_0 * exp(-(h-h0)/H)
    h0_vect = [0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 ...
        180 200 250 300 350 400 450 500 600 700 800 900 1000]';
    rho0_vect = [1.225 3.899*1e-2 1.774*1e-2 3.972*1e-3 1.057*1e-3 ...
        3.206*1e-4 8.770*1e-5 1.905*1e-5 3.396*1e-6 5.297*1e-7 9.661*1e-8 ...
        2.438*1e-8 8.484*1e-9 3.845*1e-9 2.070*1e-9 5.464*1e-10 ...
        2.789*1e-10 7.248*1e-11 2.418*1e-11 9.158*1e-12 3.725*1e-12 ...
        1.585*1e-12 6.967*1e-13 1.454*1e-13 3.614*1e-14 1.170*1e-14 ...
        5.245*1e-15 3.019*1e-15]';
    H_vect = [7.249 6.349 6.682 7.554 8.382 7.714 6.549 5.799 5.382 5.877 ...
        7.263 9.473 12.636 16.149 22.523 29.740 37.105 45.546 53.628 ...
        53.298 58.515 60.828 63.822 71.835 88.667 124.64 181.05 268.00]';

    %% Magnetic perturbation

    % Residual dipole of the spacecraft [A*m^2]
    dip_sc = [0.05 0.05 0.05];

    g01 = -29404.8;
    g11 = -1450.9;
    h11 = 4652.5;
    % H0
    H0 = sqrt(g01^2 + g11^2 + h11^2);

    %% Output
    pert_data.Fe = Fe;
    pert_data.h0_vect = h0_vect;
    pert_data.rho0_vect = rho0_vect;
    pert_data.H_vect = H_vect;
    pert_data.dip_sc = dip_sc;
    pert_data.H0 = H0;

end