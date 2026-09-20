function plot_stability_maps(alphaValues,betaValues,log10SpectralRadius, ...
    L,g,zeta,Omega,phiName,maximumFrequencyRatioShown, ...
    outputFolder,saveFigures)
%PLOT_STABILITY_MAPS Plotting and figure-export code only.

colorLimits = [-1 1];
colorMap = turbo(256);

systemDescription = sprintf([ ...
    'L=%.4g m, g=%.4g m/s^2, \\Omega=%.4g rad/s, \\zeta=%.4g, ' ...
    'u(t)=a_0\\phi(\\omega t), a_0=A\\omega^2, \\phi(\\tau)=%s'], ...
    L,g,Omega,zeta,phiName);

%% Canonical alpha-beta map

canonicalFigure = figure('Color','w');
canonicalAxes = axes(canonicalFigure);

imagesc(canonicalAxes,alphaValues,betaValues,log10SpectralRadius);
set(canonicalAxes,'YDir','normal');
hold(canonicalAxes,'on');

colormap(canonicalAxes,colorMap);
clim(canonicalAxes,colorLimits);
canonicalColorbar = colorbar(canonicalAxes);
canonicalColorbar.Label.String = ...
    'log_{10} spectral radius, log_{10} \rho(M)';

contour(canonicalAxes,alphaValues,betaValues,log10SpectralRadius, ...
    [0 0],'k','LineWidth',1.8);

xlabel(canonicalAxes,'\alpha = -\Omega^2/\omega^2');
ylabel(canonicalAxes,'\beta = a_0/(L\omega^2)');
title(canonicalAxes,{'Damped Kapitza Floquet map',systemDescription}, ...
    'Interpreter','tex');

grid(canonicalAxes,'on');
box(canonicalAxes,'on');
axis(canonicalAxes,'tight');
hold(canonicalAxes,'off');
setReadableAspect(canonicalFigure,canonicalAxes,alphaValues,betaValues);

%% Dimensional frequency-amplitude map

physicalColumns = find(alphaValues < 0);
frequencyRatio = 1./sqrt(-alphaValues(physicalColumns));

if isfinite(maximumFrequencyRatioShown)
    columnsToShow = frequencyRatio <= maximumFrequencyRatioShown;
    physicalColumns = physicalColumns(columnsToShow);
    frequencyRatio = frequencyRatio(columnsToShow);
end

physicalStabilityMap = log10SpectralRadius(:,physicalColumns);
frequencyHz = frequencyRatio*Omega/(2*pi);

% For sinusoidal support motion, beta = A/L, where A is the base
% displacement amplitude and a0 = A*omega^2 is acceleration amplitude.
amplitudeMeters = betaValues*L;
[amplitudeScale,amplitudeUnit] = readableLengthScale(amplitudeMeters);
displayedAmplitude = amplitudeMeters*amplitudeScale;

[frequencyMesh,amplitudeMesh] = meshgrid(frequencyHz,displayedAmplitude);

physicalFigure = figure('Color','w');
physicalAxes = axes(physicalFigure);

surface(physicalAxes,frequencyMesh,amplitudeMesh, ...
    zeros(size(physicalStabilityMap)),physicalStabilityMap, ...
    'EdgeColor','none');
view(physicalAxes,2);
hold(physicalAxes,'on');

colormap(physicalAxes,colorMap);
clim(physicalAxes,colorLimits);
physicalColorbar = colorbar(physicalAxes);
physicalColorbar.Label.String = ...
    'log_{10} spectral radius, log_{10} \rho(M)';

contour(physicalAxes,frequencyMesh,amplitudeMesh,physicalStabilityMap, ...
    [0 0],'k','LineWidth',1.8);

xlabel(physicalAxes,'excitation frequency, f [Hz]');
ylabel(physicalAxes,sprintf( ...
    'base displacement amplitude, A [%s]',amplitudeUnit));
title(physicalAxes,{'Physical Kapitza frequency-amplitude map', ...
    systemDescription},'Interpreter','tex');

grid(physicalAxes,'on');
box(physicalAxes,'on');
axis(physicalAxes,'tight');
hold(physicalAxes,'off');
setReadableAspect(physicalFigure,physicalAxes, ...
    frequencyHz,displayedAmplitude);

%% Save figures

if saveFigures
    exportgraphics(canonicalFigure,fullfile(outputFolder, ...
        'kapitza_stability_map_alpha_beta.png'),'Resolution',200);
    exportgraphics(physicalFigure,fullfile(outputFolder, ...
        'kapitza_stability_map_frequency_amplitude.png'),'Resolution',200);
end
end

function setReadableAspect(figureHandle,axesHandle,xValues,yValues)
% Select a useful wide aspect ratio without extreme visual distortion.
xSpan = max(xValues)-min(xValues);
ySpan = max(yValues)-min(yValues);
rawRatio = xSpan/max(ySpan,eps);
targetRatio = min(max(rawRatio,1.45),1.90);

height = 650;
width = round(height*targetRatio + 150); % Extra space for the colorbar.
position = figureHandle.Position;
position(3:4) = [width height];
figureHandle.Position = position;
pbaspect(axesHandle,[targetRatio 1 1]);
end

function [scale,unit] = readableLengthScale(valuesInMeters)
% Select a readable dimensional unit for the displacement axis.
largestValue = max(abs(valuesInMeters),[],'all');

if largestValue >= 1
    scale = 1;
    unit = 'm';
elseif largestValue >= 1e-2
    scale = 1e2;
    unit = 'cm';
elseif largestValue >= 1e-5
    scale = 1e3;
    unit = 'mm';
else
    scale = 1e6;
    unit = '\mum';
end
end
