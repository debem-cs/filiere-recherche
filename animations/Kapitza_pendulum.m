%% ====================================================================
%  KAPITZA PENDULUM  -  simulation + animation
%
%  theta measured from the UPWARD vertical.
%  Pivot:  y_p(t) = A*cos(w*t)
%  EOM:    theta'' + gamma*theta' - ((g - A*w^2*cos(w*t))/L)*sin(theta) = 0
%  Criterion: (A*w)^2 > 2*g*L
% ====================================================================
clear; close all; clc;

%% ------------------------- Parameters -------------------------------
g     = 9.81;      % gravity                              [m/s^2]
L     = 0.25;      % pendulum length                     [m]
zeta = 0.01;
gamma = 2*zeta*sqrt(g/L);
A     = 0.07;     % pivot amplitude                      [m]
f_exc = 6.0;        % pivot frequency                      [Hz]

theta0 = 10*pi/180;     % initial angle from upward vertical   [rad]
omega0 = 0.0;      % initial angular velocity             [rad/s]

t_end  = 100.0;      % total simulated time                 [s]
t_anim = 40.0;      % animated portion                     [s]
inicio_anim = 10;
fim_anim = 60;

pts_per_period = 8;     % time samples per excitation period
fps            = 60;    % playback rate                   [1/s]
trace_len      = 80;    % bob trail length                [samples]

w  = 2*pi*f_exc;
T  = 2*pi/w;
dt = 1/fps;

%% ---------------------- Stability criterion -------------------------
ratio = (A*w)^2 / (2*g*L);
fprintf('(A*w)^2 = %.3f   2*g*L = %.3f   ratio = %.2f\n', ...
        (A*w)^2, 2*g*L, ratio);
if ratio > 1
    fprintf('Critério aproximado de alta frequência satisfeito.\n');
else
    fprintf('Critério aproximado de alta frequência não satisfeito.\n');
end

%% ---------------------------- Solve ---------------------------------
tspan  = 0:dt:t_end;
odefun = @(t,x)[ x(2);
                -gamma*x(2) + ((g - A*w^2*cos(w*t))/L)*sin(x(1)) ];

% MaxStep is essential: without it the solver steps over whole
% excitation periods and reports spurious stability.
opts = odeset( ...
    'RelTol', 1e-9, ...
    'AbsTol', 1e-11, ...
    'MaxStep', T/100);
[t,x] = ode45(odefun, tspan, [theta0; omega0], opts);

theta = x(:,1);
yp    = A*cos(w*t);
xb    = L*sin(theta);
yb    = yp + L*cos(theta);

indices_anim = find(t >= inicio_anim & t <= fim_anim);

%% ------------------------- Figure set-up ----------------------------
fig = figure('Color','w','Position',[100 100 1000 460]);

% --- left panel: animation ---
ax1 = subplot(1,2,1);
lim = 1.35*L;
hold on; box on; grid on; axis equal;
set(ax1, ...
    'XLim', [-1.2*L, 1.2*L], ...
    'YLim', [-1.2*(L+A), 1.2*(L+A)]);
xlabel('x [m]'); ylabel('y [m]');

plot([0 0],[-A A],'-','Color',[0.85 0.85 0.85],'LineWidth',5);
plot([-lim lim],[0 0],'k-','LineWidth',0.5);

h_trace = plot(NaN,NaN,'-','Color',[0.90 0.60 0.25],'LineWidth',1);
h_rod   = plot(NaN,NaN,'-','Color',[0.15 0.15 0.15],'LineWidth',2.5);
h_bob   = plot(NaN,NaN,'o','MarkerSize',13,'MarkerFaceColor',[0.80 0.25 0.20], ...
               'MarkerEdgeColor','k');
h_piv   = plot(NaN,NaN,'s','MarkerSize',11,'MarkerFaceColor',[0.25 0.40 0.75], ...
               'MarkerEdgeColor','k');
h_ttl   = title('');

% --- right panel: angle history ---
ax2 = subplot(1,2,2);
hold on; box on; grid on;
plot(t, theta, '-', 'Color',[0.60 0.60 0.60],'LineWidth',0.9);
plot([0 t_end],[0 0],'k:');
set(ax2,'XLim',[0 t_end]);
xlabel('t [s]'); ylabel('theta [rad]');
title('Angle from the upward vertical');
h_mark = plot(NaN,NaN,'o','MarkerSize',7,'MarkerFaceColor',[0.80 0.25 0.20], ...
              'MarkerEdgeColor','k');

%% -------------------------- Animate ---------------------------------
fprintf('Animating %d frames...\n', numel(indices_anim));

for k = indices_anim(:).'

    if ~ishghandle(fig), break; end

    i0 = max(indices_anim(1), k-trace_len);

    set(h_trace,'XData',xb(i0:k),  'YData',yb(i0:k));
    set(h_rod,  'XData',[0 xb(k)], 'YData',[yp(k) yb(k)]);
    set(h_bob,  'XData',xb(k),     'YData',yb(k));
    set(h_piv,  'XData',0,         'YData',yp(k));
    set(h_mark, 'XData',t(k),      'YData',theta(k));
    set(h_ttl,  'String', sprintf('t = %.3f s    theta = %+.3f rad', ...
                                   t(k), theta(k)));

    drawnow;
    pause(1/fps);
end

fprintf('Done. Final angle: %+.4f rad at t = %.2f s\n', theta(end), t(end));