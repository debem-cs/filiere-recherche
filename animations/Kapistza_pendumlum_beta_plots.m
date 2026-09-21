%% ====================================================================
%  TILTED-EXCITATION PENDULUM  (Ciezkowski, JSV 491, 2021)
%
%  theta measured from the upward vertical.
%  Pivot:  r_p(t) = A*cos(w*t) * [sin(beta); cos(beta)]
%  EOM:    L*theta'' = g*sin(theta) + A*w^2*cos(w*t)*sin(beta - theta)
%                      - L*gamma*theta'
%  beta = 0 recovers the classical Kapitza pendulum.
%  Stability conditions: Eq. (6) of the paper.
% ====================================================================
clear; close all; clc;

%% Parameters
g     = 9.81;                   % gravity [m/s^2]
L     = 0.25;                   % pendulum length [m]
zeta  = 0.01;                   % damping ratio
gamma = 2*zeta*sqrt(g/L);       % damping coefficient [1/s]
A     = 0.09;                   % pivot amplitude [m]
f_exc = 6.0;                    % pivot frequency [Hz]

phi_s  = deg2rad(60);           % desired equilibrium [rad]
theta0 = phi_s + deg2rad(10);   % initial angle [rad]
omega0 = 0.0;                   % initial angular velocity [rad/s]

t_end      = 100.0;             % simulated time [s]
start_anim = 20.0;              % animation window [s]
end_anim   = 60.0;

fps       = 60;                 % sampling / playback rate [1/s]
trace_len = 30;                 % bob trail length [samples]
save_png  = true;               % export final figure

w  = 2*pi*f_exc;
T  = 2*pi/w;
dt = 1/fps;

%% Stability check and beta
lambda = (A*w)^2 / (2*g*L);     % analogue of 3*A^2*Om^2/(4*g*l) in the paper

if phi_s <= pi/2
    lambda_min = sqrt(4*sin(phi_s)^2 + cos(phi_s)^2);   % Eq. (6)
else
    lambda_min = 2*sin(phi_s);                          % Eq. (6)
end

fprintf('lambda          = %.4f\n', lambda);
fprintf('lambda required = %.4f   (phi_s = %.1f deg)\n', ...
        lambda_min, rad2deg(phi_s));

if lambda <= lambda_min
    error(['lambda = %.3f is below the required %.3f for phi_s = %.1f deg. ' ...
           'Increase A or f_exc, or move phi_s closer to 0.'], ...
          lambda, lambda_min, rad2deg(phi_s));
end

beta = 0.5*(2*phi_s - asin(2*sin(phi_s)/lambda));      % Eq. (6)
fprintf('beta            = %.4f rad = %.2f deg\n', beta, rad2deg(beta));

% Slow-motion natural frequency (curvature of the effective potential)
w0_sq = (A^2*w^2/(2*L^2))*cos(2*(beta - phi_s)) - (g/L)*cos(phi_s);
if w0_sq <= 0
    error('No potential minimum at phi_s (w0^2 = %.4f).', w0_sq);
end
w0 = sqrt(w0_sq);
fprintf('slow frequency  = %.3f rad/s  -> period %.4f s\n\n', w0, 2*pi/w0);

%% Solve
tspan  = 0:dt:t_end;
odefun = @(t,x)[ x(2);
                -gamma*x(2) + (g*sin(x(1)) ...
                               + A*w^2*cos(w*t)*sin(beta - x(1)))/L ];

% MaxStep prevents the solver from stepping over excitation periods
opts  = odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',T/100);
[t,x] = ode45(odefun, tspan, [theta0; omega0], opts);

theta = x(:,1);
xp    = A*cos(w*t)*sin(beta);   % pivot
yp    = A*cos(w*t)*cos(beta);
xb    = xp + L*sin(theta);      % bob
yb    = yp + L*cos(theta);

%% Position error
err_deg = rad2deg(phi_s - theta);

% Signed moving average over one excitation period (edge-normalized).
% Removes the fast component and leaves the slow error phi_s - phi.
N       = max(1, round(T/dt));
k_win   = ones(N,1);
err_avg = conv(err_deg, k_win, 'same') ./ conv(ones(size(err_deg)), k_win, 'same');

frames = find(t >= start_anim & t <= end_anim);
if isempty(frames)
    error('Animation window [%.1f, %.1f] s lies outside the simulation.', ...
          start_anim, end_anim);
end

%% Figure
fig = figure('Color','w','Position',[80 80 1150 600]);

% Animation
ax1 = subplot(2,2,[1 3]);
lim = 1.35*L + A;
hold on; box on; grid on; axis equal;
set(ax1,'XLim',[-lim lim],'YLim',[-lim lim]);
xlabel('x [m]'); ylabel('y [m]');

plot([-lim lim],[0 0],'k-','LineWidth',0.5);
plot([0 0],[-lim lim],'k-','LineWidth',0.5);
plot(A*[-1 1]*sin(beta), A*[-1 1]*cos(beta), ...           % pivot track
     '-','Color',[0.85 0.85 0.85],'LineWidth',5);
plot([0 1.25*L*sin(phi_s)],[0 1.25*L*cos(phi_s)], ...      % target direction
     '--','Color',[0.30 0.65 0.30],'LineWidth',1.2);

h_trace = plot(NaN,NaN,'-','Color',[0.90 0.60 0.25],'LineWidth',1);
h_rod   = plot(NaN,NaN,'-','Color',[0.15 0.15 0.15],'LineWidth',2.5);
h_bob   = plot(NaN,NaN,'o','MarkerSize',13, ...
               'MarkerFaceColor',[0.80 0.25 0.20],'MarkerEdgeColor','k');
h_piv   = plot(NaN,NaN,'s','MarkerSize',11, ...
               'MarkerFaceColor',[0.25 0.40 0.75],'MarkerEdgeColor','k');
h_ttl   = title('');

% Angle history
ax2 = subplot(2,2,2);
hold on; box on; grid on;
h_th  = plot(t, rad2deg(theta), '-', 'Color',[0.60 0.60 0.60],'LineWidth',0.9);
h_tgt = plot([0 t_end], rad2deg(phi_s)*[1 1], '--', 'Color',[0.30 0.65 0.30]);
set(ax2,'XLim',[0 t_end]);
xlabel('t [s]'); ylabel('\theta [deg]');
title(sprintf('\\beta = %.1f deg   \\rightarrow   \\phi_s = %.1f deg', ...
      rad2deg(beta), rad2deg(phi_s)));
h_mark1 = plot(NaN,NaN,'o','MarkerSize',7, ...
               'MarkerFaceColor',[0.80 0.25 0.20],'MarkerEdgeColor','k');
legend([h_th h_tgt], {'\theta(t)', 'Target \phi_s'}, 'Location','best');

% Position error
ax3 = subplot(2,2,4);
hold on; box on; grid on;
h_e    = plot(t, err_deg, '-', 'Color',[0.60 0.60 0.60],'LineWidth',0.9);
h_eavg = plot(t, err_avg, '-', 'Color',[0.25 0.40 0.75],'LineWidth',1.6);
plot([0 t_end],[0 0],'k:');
set(ax3,'XLim',[0 t_end]);
xlabel('t [s]'); ylabel('e [deg]');
title('Position error  e = \phi_s - \theta');
h_mark2 = plot(NaN,NaN,'o','MarkerSize',7, ...
               'MarkerFaceColor',[0.80 0.25 0.20],'MarkerEdgeColor','k');
legend([h_e h_eavg], {'e(t)', 'Period average'}, 'Location','northeast');

%% Animate
fprintf('Animating %d frames (t = %.1f to %.1f s)...\n', ...
        numel(frames), t(frames(1)), t(frames(end)));

for k = frames(:).'
    if ~ishghandle(fig), break; end
    i0 = max(1, k - trace_len);

    set(h_trace,'XData',xb(i0:k),      'YData',yb(i0:k));
    set(h_rod,  'XData',[xp(k) xb(k)], 'YData',[yp(k) yb(k)]);
    set(h_bob,  'XData',xb(k),         'YData',yb(k));
    set(h_piv,  'XData',xp(k),         'YData',yp(k));
    set(h_mark1,'XData',t(k),          'YData',rad2deg(theta(k)));
    set(h_mark2,'XData',t(k),          'YData',err_deg(k));
    set(h_ttl,  'String', sprintf('t = %.2f s    \\theta = %.1f deg', ...
                                   t(k), rad2deg(theta(k))));
    drawnow;
    pause(1/fps);
end

%% Diagnostics
idx = t > t_end - 10;
fprintf('\nMean angle over the last 10 s : %.2f deg   (target %.2f deg)\n', ...
        rad2deg(mean(theta(idx))), rad2deg(phi_s));
fprintf('Std  over the same window     : %.2f deg\n', ...
        rad2deg(std(theta(idx))));
fprintf('Mean error over the same window: %.3f deg\n', mean(err_deg(idx)));

if save_png && ishghandle(fig)
    print(fig, '-dpng', '-r300', 'tilted_pendulum.png');
    fprintf('Figure saved to tilted_pendulum.png\n');
end