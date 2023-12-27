function sc_data = spaceCraftData()
    % Output a struct containing all the data relative to the spacecraft:
    %   - I_mat:    Principal moments of Inertia matrix [Kg*m^2]
    %   - I_inv:    Inverse of principal inertia matrix [(Kg*m^2)^-1]
    %   - Ir:       Inertia moment of inertia wheel [Kg*m^2]
    %   - NB:       Normals to the body [-]
    %   - rhoS:     Specular reflection coefficient [-]
    %   - rhoD:     Diffusione reflection coefficient [-]
    %   - A:        Surface area [m^2]
    %   - rF:       Distance of body parts from baricentre [m]

    % Principal moments of Intertia [Kg*m^2]
    % I = [100.9 25.1 91.6] * 1e-2;
    % I = [0.04 0.06 0.08];
    % I = [279 945 1085];
    I = [4.2929 13.2152 13.6491];
    % I = [0.0700 0.0504 0.0109];
    % Principal inertia matrix [Kg*m^2]
    I_mat = I.*eye(3);
    % Inverse of intertia matrix
    I_inv = inv(I_mat);
    % NOTE: The inverse is calculated before hand to save computation time
    % as the inertia matrix does not change during the mission

    % Inertia moment of inertia wheel [Kg*m^2]
    Ir = [0 3.183e-2 0];

    % Normals to body 
    NB = [ 1  0  0;      % Body 1       
           0  1  0;      % Body 2
          -1  0  0;      % Body 3
           0 -1  0;      % Body 4
           0  0  1;      % Body 5
           0  0 -1;      % Body 6
           1  0  0;      % Panel 1
          -1  0  0;      % Panel 2
           1  0  0;      % Panel 3
          -1  0  0];     % Panel 4

    % specular refelction coefficient
    rhoS_body = 0.5;
    rhoS_panel = 0.1;
    rhoS = [repmat(rhoS_body, 6, 1); repmat(rhoS_panel, 4, 1)];
    
    % diffusion reflection coefficient
    rhoD = repmat(0.1, 10, 1);

    % Surface area [m^2]
    % A_body = [6e-2*ones(4, 1); 4e-2*ones(2, 1)];
    A_body = [0.64;0.4;0.64;0.4; 0.4*ones(2,1)];
    A_panel = 12e-2*ones(4,1);
    A_sc = [A_body; A_panel];

    % distance from baricentre [cm]
    rF = [ 10    0   0;       % Body 1
            0   10   0;       % Body 2
          -10    0   0;       % Body 3
            0  -10   0;       % Body 4
            0    0   15;      % Body 5
            0    0  -15;      % Body 6
            0    0   45;      % Panel 1
            0    0   45;      % Panel 2
            0    0  -45;      % Panel 3
            0    0  -45];     % Panel 4
    % convert to meters
    rF = rF * 1e-2;

    % Drag coefficient
    cD = 2.1;

    sc_data.I_mat = I_mat;
    sc_data.I_inv = I_inv;
    sc_data.Ir = Ir;
    sc_data.NB = NB;
    sc_data.rhoS = rhoS;
    sc_data.rhoD = rhoD;
    sc_data.A = A_sc;
    sc_data.rF = rF;
    sc_data.magnetorquer = 80; %[Am^2] linear dipole moment 
    sc_data.wr = (3600*60)/6.28; %[rad/s] nominal rotation speed
    sc_data.NumFaces = size(NB, 1);
    sc_data.cD = cD;
end