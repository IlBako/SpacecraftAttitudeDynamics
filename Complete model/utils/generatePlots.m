if exist("de_tumb", 'var')
    out = de_tumb;
elseif exist("no_cont", 'var')
    out = no_cont;
elseif exist("point", 'var')
    out = point;
end

% Plot - om
figure()
plot(out.time, out.dynamics_omega, 'LineWidth', 1);
grid on;
title("Angular velocities");
xlabel("Time\ [s]", 'Interpreter','latex');
ylabel("$\omega$ [rad/s]", 'Interpreter','latex');
legend("$\omega_x$", "$\omega_y$", "$\omega_z$", 'Location','northeast', 'Interpreter', 'latex');

if alg_idx == 2
figure()
plot(out.time, out.Mc, 'LineWidth', 1);
title("Control moment");
xlabel("Time \[s]", 'Interpreter', 'latex');
ylabel("$M_c\ [N*m]$", 'Interpreter','latex');

figure();
subplot(2,1,1);
plot(out.time, out.Mc);%, 'LineWidth', 1);
title("Control moment");
xlabel("Time \[s]", 'Interpreter', 'latex');
ylabel("$M_c\ [N*m]$", 'Interpreter','latex');
legend("$M_{c\_x}$", "$M_{c\_y}$", "$M_{c\_z}$", 'Interpreter', 'latex');

subplot(2,1,2);
plot(out.time, out.MC_actuators);
title("Actuator control moment");
xlabel("Time \[s]", 'Interpreter', 'latex');
ylabel("$M_c\ [N*m]$", 'Interpreter','latex');
legend("$M_c\_x$", "$M_c\_y$", "$M_c\_z$", 'Interpreter', 'latex');
end

% Plot - disturbs
figure()
plot(out.time, vecnorm(out.T_GG, 2, 2), 'LineWidth', 1, 'DisplayName', 'Gravity Gradient');
grid on, hold on
plot(out.time, vecnorm(out.T_Mag, 2, 2), 'LineWidth', 1, 'DisplayName', 'Magnetic');
plot(out.time, vecnorm(out.T_aero, 2, 2), 'LineWidth', 1, 'DisplayName', 'Aerodynamic drag');
xlabel("Time\ [s]", 'Interpreter','latex');
ylabel("T $[N*m]$", 'Interpreter','latex');
legend('Location','northeast', 'Interpreter', 'latex');
title("Disturbances moment");

figure()
subplot(3,1,1)
plot(out.time, vecnorm(out.T_GG, 2, 2), 'LineWidth', 1);
title("Gravity Gradient");
subplot(3,1,2)
plot(out.time, vecnorm(out.T_Mag, 2, 2), 'LineWidth', 1);
title("Magnetic");
subplot(3,1,3)
plot(out.time, vecnorm(out.T_aero, 2, 2), 'LineWidth', 1);
title("Aerodynamic drag");

% Plot - kinetic energy
figure()
T = 0.5 * (sc_data.I_mat(1,1)*out.dynamics_omega(:,1).^2 + sc_data.I_mat(2,2)*out.dynamics_omega(:,2).^2 + sc_data.I_mat(3,3)*out.dynamics_omega(:,3).^2);
plot(out.time, T, 'LineWidth', 1)
grid on
title("Kinetic energy")
xlabel("Time\ [s]", 'Interpreter','latex');
ylabel("T [J]", 'Interpreter','latex');
legend("Kinetic Energy", 'Location','northeast', 'Interpreter', 'latex');

if alg_idx == 2
figure();
plot(out.time, out.dynamics_omega, 'LineWidth', 1); 
hold on; grid on; 
plot(out.time, out.kepler_th_dot, 'LineWidth', 1.5); 
for i=0:floor(out.time(end)/orbit_data.T)
    xline(out.time(1) + i*orbit_data.T, 'k--', 'LineWidth', 0.25, 'HandleVisibility','off');
end
xlabel("Time\ [s]", 'Interpreter','latex');
ylabel("$\omega\ [rad/s]$", 'Interpreter','latex');
legend('$\omega_x$', '$\omega_y$', '$\omega_z$', '$\dot{\theta}$', 'Interpreter', 'latex');
end

% Plot - err om_BL
% figure
% plot(out.time, out.w_BL - out.dynamics_omega, 'LineWidth', 1) % TODO check if the formula is correct
% grid on
% title("Angular velocity error")
% xlabel("Time\ [s]", 'Interpreter','latex');
% ylabel("$\omega$ [rad/s]", 'Interpreter','latex');
% legend("$\omega_{x} error$", "$\omega_{y} error$", "$\omega_{z} error$", 'Location','northeast', 'Interpreter', 'latex');

% Plot - stability graph
% I_x = sc_data.I_mat(1,1);
% I_y = sc_data.I_mat(2,2);
% I_z = sc_data.I_mat(3,3);
% 
% K_x = (I_z - I_y)/I_x;
% K_y = (I_z - I_x)/I_y;
% K_z = (I_y - I_x)/I_z;
% 
% figure()
% plot(K_y, K_x, 'or')
% hold on
% plot([-1 1], [-1 1])
% plot([-1 1], [0 0], color = '[0.2, 0.2, 0.2]')
% plot([0 0], [-1 1], color = '[0.2, 0.2, 0.2]')
% grid on, axis equal
% xlim([-1 1])
% ylim([-1 1])
% title("Stability plot")

% Plot - orbit
figure()
earthPlot();
plot3(out.kepler_r_vec(:,1), out.kepler_r_vec(:,2), out.kepler_r_vec(:,3), 'LineWidth', 1);
grid on, axis equal
view(30, 20)
title("Orbit")