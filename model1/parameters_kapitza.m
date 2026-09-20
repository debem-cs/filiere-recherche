%% simulation
g = 9.81;
l = 0.25;
m = 1.0;
b = 0.0628;

%change amplitude and frequency
A = 0.04;
f = 15;

%change initial condition
x0 = [75*pi/180; 0];

p = [g; l; m; b; A; f];

model_name = 'kapitza_manual';

load_system(model_name);

fprintf('Iniciando simulação de %s...\n', model_name);
out = sim(model_name);
fprintf('Simulação concluída!\n');

%% Graphic visualization
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
subtitle(sprintf('f = %.2f Hz, A = %.3f m, L = %.3f m, b = %.3f, theta_0 = %.2f°', f, A, l, b, x0(1)*180/pi))

nexttile
plot(t_omega, omega)
grid on
xlabel('Time [s]')
ylabel('Angular speed [rad/s]')