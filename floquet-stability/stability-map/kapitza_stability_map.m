%% FLOQUET STABILITY MAP OF A DAMPED KAPITZA PENDULUM
% Engineering calculation script.
%
% The suspension acceleration is
%       u(t) = a0*phi(omega*t).
%
% The nonlinear equation around the upward position theta = 0 is
%       theta_ddot + 2*zeta*Omega*theta_dot
%       + [-Omega^2 + u(t)/L]*sin(theta) = 0.
%
% After tau = omega*t, its linearization is
%       theta'' + 2*delta*theta'
%       + [alpha + beta*phi(tau)]*theta = 0,
%
% with
%       alpha = -Omega^2/omega^2,
%       beta  = a0/(L*omega^2),
%       delta = zeta*Omega/omega.

clearvars;
close all;

%% 1. PHYSICAL PENDULUM

L = 1.0;                    % Pendulum length [m]
g = 9.81;                   % Gravitational acceleration [m/s^2]
zeta = 0.01;                % Physical damping ratio [-]
Omega = sqrt(g/L);          % Natural angular frequency [rad/s]

%% 2. PERIODIC INPUT AND NONLINEAR EQUATION

syms tau theta thetaDot alpha beta delta real

% Change this expression to study another 2*pi-periodic waveform.
phi = cos(tau);
phiName = 'cos(\tau)';

% Nonlinear normalized state equation x' = f(tau,x).
x = [theta; thetaDot];
f = [thetaDot; ...
    -2*delta*thetaDot - (alpha + beta*phi)*sin(theta)];

% Linearize symbolically around the upward equilibrium theta = thetaDot = 0.
A_symbolic = simplify(subs(jacobian(f,x),x,[0;0]));

fprintf('phi(tau) = %s\n',char(phi));
fprintf('Linearized state matrix A(tau):\n');
disp(A_symbolic);

%% 3. STABILITY-MAP DOMAIN

% These are the quantities selected by the engineer.
alphaValues = linspace(-1.0,0.25,1000);
betaValues = linspace(0.0,4.0,300);

% Number of constant-coefficient intervals used during one forcing period.
% Doubling this value provides a simple time-discretization convergence test.
numberOfTimeSteps = 120;

% Only affects the dimensional frequency-amplitude figure. Set to Inf to
% display every frequency represented by the negative-alpha grid points.
maximumFrequencyRatioShown = 50.0;   % Maximum omega/Omega shown

showProgress = true;
saveResults = true;
saveFigures = true;

%% 4. DISCRETIZE ONE PERIOD OF THE INPUT

period = 2*pi;
timeStep = period/numberOfTimeSteps;
timeMidpoints = ((0:numberOfTimeSteps-1) + 0.5)*timeStep;

% The waveform does not depend on alpha or beta, so evaluate it only once.
phiAtMidpoints = double(subs(phi,tau,timeMidpoints));

%% 5. CALCULATE THE MONODROMY MATRIX AT EVERY GRID POINT

numberOfAlphaValues = numel(alphaValues);
numberOfBetaValues = numel(betaValues);

log10SpectralRadius = nan(numberOfBetaValues,numberOfAlphaValues);
determinantRelativeError = nan(numberOfBetaValues,numberOfAlphaValues);

numberOfGridPoints = numberOfAlphaValues*numberOfBetaValues;
completedGridPoints = 0;
nextProgressReport = 0.1;
calculationTimer = tic;

for alphaIndex = 1:numberOfAlphaValues
    alphaValue = alphaValues(alphaIndex);

    % For alpha < 0, delta = zeta*sqrt(-alpha) follows directly from the
    % physical definitions. abs(alpha) extends the same damping convention
    % to the alpha > 0 portion of the canonical diagram.
    deltaValue = zeta*sqrt(abs(alphaValue));

    % Liouville's formula gives det(M) exactly for this system.
    expectedDeterminant = exp(-2*deltaValue*period);

    for betaIndex = 1:numberOfBetaValues
        betaValue = betaValues(betaIndex);

        monodromyMatrix = calculateMonodromyMatrix( ...
            alphaValue,betaValue,deltaValue, ...
            phiAtMidpoints,timeStep);

        floquetMultipliers = eig(monodromyMatrix);
        spectralRadius = max(abs(floquetMultipliers));

        log10SpectralRadius(betaIndex,alphaIndex) = ...
            log10(max(spectralRadius,realmin));

        determinantRelativeError(betaIndex,alphaIndex) = ...
            abs(det(monodromyMatrix)-expectedDeterminant) ...
            / max(expectedDeterminant,realmin);

        completedGridPoints = completedGridPoints + 1;
        completedFraction = completedGridPoints/numberOfGridPoints;

        if showProgress && completedFraction >= nextProgressReport
            fprintf('%3.0f %% complete, elapsed %.1f s\n', ...
                100*completedFraction,toc(calculationTimer));
            nextProgressReport = nextProgressReport + 0.1;
        end
    end
end

stable = log10SpectralRadius < 0;
maximumDeterminantError = max(determinantRelativeError,[],'all');

fprintf('Maximum Liouville determinant relative error: %.3e\n', ...
    maximumDeterminantError);

%% 6. PLOT AND SAVE THE RESULTS

projectFolder = fileparts(mfilename('fullpath'));
outputFolder = fullfile(projectFolder,'generated');

if (saveResults || saveFigures) && ~isfolder(outputFolder)
    mkdir(outputFolder);
end

plot_stability_maps(alphaValues,betaValues,log10SpectralRadius, ...
    L,g,zeta,Omega,phiName,maximumFrequencyRatioShown, ...
    outputFolder,saveFigures);

if saveResults
    save(fullfile(outputFolder,'kapitza_stability_data.mat'), ...
        'alphaValues','betaValues','log10SpectralRadius','stable', ...
        'determinantRelativeError','maximumDeterminantError', ...
        'L','g','zeta','Omega','numberOfTimeSteps','phiName');
end

%% LOCAL NUMERICAL FUNCTION

function M = calculateMonodromyMatrix(alpha,beta,delta,phiValues,h)
%CALCULATEMONODROMYMATRIX Explicit transition product over one period.
%
% At each midpoint the coefficient matrix is treated as constant:
%
%       A_k = [ 0,                         1
%              -(alpha + beta*phi_k),  -2*delta ].
%
% The exact transition of that frozen system is E_k = exp(A_k*h). The
% transition matrices are multiplied in chronological order:
%
%       M = E_(N-1)*...*E_1*E_0.
%
% A closed 2-by-2 formula is used instead of MATLAB's general expm function.

q = alpha + beta*phiValues;
discriminant = delta^2-q;

cosineTerm = zeros(size(discriminant));
sineTerm = zeros(size(discriminant));

hyperbolicSteps = discriminant > 1e-13;
oscillatorySteps = discriminant < -1e-13;
repeatedRootSteps = ~(hyperbolicSteps | oscillatorySteps);

root = sqrt(discriminant(hyperbolicSteps));
cosineTerm(hyperbolicSteps) = cosh(root*h);
sineTerm(hyperbolicSteps) = sinh(root*h)./root;

root = sqrt(-discriminant(oscillatorySteps));
cosineTerm(oscillatorySteps) = cos(root*h);
sineTerm(oscillatorySteps) = sin(root*h)./root;

% lim sinh(root*h)/root = h when root tends to zero.
cosineTerm(repeatedRootSteps) = 1;
sineTerm(repeatedRootSteps) = h;

decay = exp(-delta*h);

% Start with the identity transition at tau = 0.
m11 = 1;
m12 = 0;
m21 = 0;
m22 = 1;

for timeIndex = 1:numel(phiValues)
    % Exact exponential of the frozen 2-by-2 state matrix.
    e11 = decay*(cosineTerm(timeIndex) ...
        + delta*sineTerm(timeIndex));
    e12 = decay*sineTerm(timeIndex);
    e21 = -decay*q(timeIndex)*sineTerm(timeIndex);
    e22 = decay*(cosineTerm(timeIndex) ...
        - delta*sineTerm(timeIndex));

    % Left multiplication preserves chronological order.
    next11 = e11*m11 + e12*m21;
    next12 = e11*m12 + e12*m22;
    next21 = e21*m11 + e22*m21;
    next22 = e21*m12 + e22*m22;

    m11 = next11;
    m12 = next12;
    m21 = next21;
    m22 = next22;
end

M = [m11 m12; m21 m22];
end
