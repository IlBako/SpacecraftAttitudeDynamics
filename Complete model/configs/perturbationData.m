function pert_data = perturbationData(n)  
    % Output a struct containing all the data relative to the perturbations:
    % GG Perturbation:
    %   - w_LN:     Reference moving angular velocity [rad/s]
    % SRP Perturbation
    %   - Fe:       Solar radiation intensity [W/m^2]
    % Magnetic Perturbation:
    %   - H0:        
    %   - m:
    % Required inputs:
    %   -n:         

    % Reference moving angular velocity [rad/s]
    w_LN = [0; 0; n];                   

    % Solar radiation intensity [W/m^2]
    Fe = 1358;

    g01 = -29404.8;
    g11 = -1450.9;
    h11 = 4652.5;
    % H0
    H0 = sqrt(g01^2 + g11^2 + h11^2);

    pert_data.w_LN = w_LN;
    pert_data.Fe = Fe;
    pert_data.H0 = H0;

end