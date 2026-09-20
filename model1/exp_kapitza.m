parameters_kapitza

modelo = 'kapitza_manual';

% choose different amplitudes and frequencies to test
casos = [
    0.04  7
    0.04  15
    0.02  15
    0.04  20
];

resultados = cell(size(casos,1),1);

figure
tiledlayout(2,2)

for k = 1:size(casos,1)

    Ak = casos(k,1);
    fk = casos(k,2);

    pk = [g; l; m; b; Ak; fk];

    entrada = Simulink.SimulationInput(modelo);

    entrada = entrada.setVariable('p', pk);
    entrada = entrada.setVariable('x0', x0);

    entrada = entrada.setModelParameter( ...
        'StopTime', '100', ...
        'MaxStep', num2str(1/(100*fk),17));

    saida = sim(entrada);

    t_theta = saida.theta_ts.Time(:);
    theta = saida.theta_ts.Data(:);

    t_omega = saida.omega_ts.Time(:);
    omega = saida.omega_ts.Data(:);

    assert(numel(t_theta) == numel(theta), ...
        'O tempo e o angulo têm tamanhos diferentes.');

    assert(numel(t_omega) == numel(omega), ...
        'O tempo e a velocidade têm tamanhos diferentes.');

    resultados{k} = struct( ...
        'A', Ak, ...
        'f', fk, ...
        't_theta', t_theta, ...
        'theta', theta, ...
        't_omega', t_omega, ...
        'omega', omega);

    nexttile
    plot(t_theta, theta)
    yline(0, ':')
    grid on

    xlabel('Time [s]')
    ylabel('\theta [degrees]')
    title(sprintf('A = %.3f m; f = %.1f Hz', Ak, fk))
end