clear; close all; clc;

% Physical parameters
g = 9.81;
L = 0.25;

lambda = 1;
zeta = lambda/(2*sqrt(g/L));

% Vibrational control parameters
w = 500;
f = w/(2*pi);

alpha_ctrl = 40;
A = 2*alpha_ctrl/w^(3/2);

% Controller gains
c1 = 1;
c2 = 1.5;

% The LgV controller requires the third element
pc = [c1; c2; alpha_ctrl];
p = [g; L; zeta; A; f];

% Same initial condition for both cases
x0 = [150*pi/180; 0];

% % Models
modelos = {'model2_control', 'model2_control_LgV'};

% Same numerical settings
for k = 1:numel(modelos)

    load_system(modelos{k});

    set_param(modelos{k}, ...
        'SolverType', 'Variable-step', ...
        'Solver', 'ode45', ...
        'RelTol', '1e-7', ...
        'AbsTol', '1e-9', ...
        'MaxStep', num2str(1/(100*f),17), ...
        'StopTime', '10');
end

% Store the results separately
fprintf('Simulating the vibrational controller...\n');
out_vib = sim(modelos{1});

fprintf('Simulating the non-vibrational LgV controller...\n');
out_lgv = sim(modelos{2});

fprintf('Simulations completed.\n');

% Angles: already in degrees in both models
theta_vib = out_vib.theta_ts;
theta_lgv = out_lgv.theta_ts;

% Angular velocities: rad/s
vel_vib = out_vib.omega_ts;
vel_lgv = out_lgv.omega_ts;
Vx = out_vib.Vx;

% Comparison
figure;
tiledlayout(2,1);

% Angle
ax1 = nexttile;

plot(theta_vib.Time(:), theta_vib.Data(:), ...
    'b-', 'LineWidth', 1.3);
hold on;

plot(theta_lgv.Time(:), theta_lgv.Data(:), ...
    'r--', 'LineWidth', 1.3);

yline(0, 'k:', 'HandleVisibility', 'off');

grid on;
xlabel('Time [s]');
ylabel('\theta [degrees]');
title('Angle from the Upper Vertical');
legend('Vibrational', 'Non-vibrational LgV', ...
    'Location', 'best');

% Angular velocity
ax2 = nexttile;

plot(vel_vib.Time(:), vel_vib.Data(:), ...
    'b-', 'LineWidth', 1.3);
hold on;

plot(vel_lgv.Time(:), vel_lgv.Data(:), ...
    'r--', 'LineWidth', 1.3);

yline(0, 'k:', 'HandleVisibility', 'off');

grid on;
xlabel('Time [s]');
ylabel('Angular velocity [rad/s]');
title('Angular Velocity');
legend('Vibrational', 'Non-vibrational LgV', ...
    'Location', 'best');

linkaxes([ax1, ax2], 'x');

sgtitle('Controller Comparison');

% Angle comparison
figure;

plot(theta_vib.Time(:), theta_vib.Data(:), ...
    'b-', 'LineWidth', 1.3);
hold on;

plot(theta_lgv.Time(:), theta_lgv.Data(:), ...
    'r--', 'LineWidth', 1.3);

yline(0, 'k:', 'HandleVisibility', 'off');

grid on;
xlabel('Time [s]');
ylabel('\theta [degrees]');
title('Controller Comparison — Angle from the Upper Vertical');

legend('Vibrational', 'Non-vibrational LgV','Location', 'best');

% % plot3(theta_vib.Data(:), vel_vib.Data(:), Vx.Data(:), 'LineWidth', 2);
% % grid on;
% % xlabel('x1'); ylabel('x2'); zlabel('V');
% 
