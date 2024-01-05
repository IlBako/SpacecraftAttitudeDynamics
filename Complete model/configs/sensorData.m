function sensor_data = sensorData()

    %% Magnetometer 

    % [T] Full scale measurement of the sensor (range is +-)
    magnetometer.FSS = 100e-6;
    % [%] Accuracy of the sensor
    magnetometer.acc = 0.05;
    % [%] Non-linearity of the sensor
    magnetometer.lin = 0.015;
    % [V/T] Sensitivity of the sensor
    magnetometer.sens = 100*1e3;
    % [T] Standard deviation of the sensor
    magnetometer.std_dev = sqrt(magnetometer.acc ^ 2 + magnetometer.lin^2)/sqrt(3) * magnetometer.FSS / 100;
    % [T^2] Variance of the sensor
    magnetometer.variance = magnetometer.std_dev^2;
    % [V] Quantization interval in Volt
    V_quant = 0.025*2;
    % [T] Quantization in Tesla
    magnetometer.T_quant = V_quant/magnetometer.sens;
    % [Hz] Frequency of the sensor
    magnetometer.freq = 10;

    % To calculate the magnetometer non orthogonality error 
    % we use the inverse approach of hard and soft iron calibration.
    % The way hard and soft iron calibration work is by removing the offset
    % of the center of the magnetometer sphere and by transforming the
    % noisy magnetometer measure, which creates an ellipse, in a sphere

    a = deg2rad(rand(3,1)*1); % Soft iron - orthogonality within 1 deg
    b = deg2rad(rand(3,1)*360); % Hard iron - rotation over 360 deg
    % Orthogonality error matrix: this matrix both rotates and stretches
    % the magnetic field to turn it into an ellipse with offset center
    magnetometer.A_ortho = [cos(a(1))               sin(a(1))*cos(b(1))     sin(a(1))*sin(b(1));
                            sin(a(2))*cos(b(2))     cos(a(2))               sin(a(2))*sin(b(2));
                            sin(a(3))*cos(b(3))     sin(a(3))*sin(b(3))     cos(a(3))];

    %% Earth horizon sensor (range is +-)

    % [rad] Accuracy of the sensor (0.05deg at 3 sigma)
    horizon.acc = deg2rad(0.05)/3;
    % Misalignment error
    horizon.misalign = deg2rad(0.5);
    % [Hz] Frequency of the sensor
    horizon.freq = 10;
    
    %% Output

    sensor_data.magnetometer = magnetometer;
    sensor_data.horizon = horizon;

end