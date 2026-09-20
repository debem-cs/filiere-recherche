clear; close all; clc;

g = 9.81;
L = 0.25;
% zeta = 0.01;
lambda = 1;
zeta = lambda/(2*sqrt(g/L));

A = 0.08;
f = 15;

x0 = [45*pi/180; 0];

p = [g; L; zeta; A; f];

model_name = 'kapitza_manual_model2';
load_system(model_name);

% Ajustes para resolver a excitação periódica
set_param(model_name, ...
    'SolverType', 'Variable-step', ...
    'Solver', 'ode45', ...
    'RelTol', '1e-7', ...
    'AbsTol', '1e-9', ...
    'MaxStep', num2str(1/(100*f),17), ...
    'StopTime', '20');

fprintf('Iniciando simulação de %s...\n', model_name);
out = sim(model_name);
fprintf('Simulação concluída!\n');

%% visualize

% Ângulo e seu tempo
t_theta = out.theta_ts.Time(:);
theta = out.theta_ts.Data(:);

% Velocidade angular e seu tempo
t_omega = out.omega_ts.Time(:);
omega = out.omega_ts.Data(:);

% Verificações
assert(numel(t_theta) == numel(theta), ...
    'O tempo e o angulo têm tamanhos diferentes.');

assert(numel(t_omega) == numel(omega), ...
    'O tempo e a velocidade têm tamanhos diferentes.');

figure
tiledlayout(2,1)

nexttile
plot(t_theta, theta)
yline(0, ':')
grid on
xlabel('Time [s]')
ylabel('\theta [degrees]')
title('Angle from the upper vertical')
subtitle(sprintf('f = %.2f Hz, A = %.3f m, L = %.3f m, \\zeta = %.3f,theta_0 = %.2f°', f, A, L, zeta, x0(1)*180/pi))

nexttile
plot(t_omega, omega)
grid on
xlabel('Time [s]')
ylabel('Angular speed [rad/s]')