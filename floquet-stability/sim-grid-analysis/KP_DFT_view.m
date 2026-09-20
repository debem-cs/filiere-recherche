%% INTERACTIVE VIEW OF THE STORED KAPITZA DFT RESULTS
%
% Run this script and click a green point in the amplitude-frequency map.
% Press Enter while the grid figure is active to finish.

clear;
clc;
close all;

%% 1. LOAD THE RESULTS

projectFolder = fileparts(mfilename('fullpath'));
outputFolder = fullfile(projectFolder,'generated');
dftDataFile = fullfile(outputFolder,'KP_DFT_results.mat');

if ~isfile(dftDataFile)
    error(['DFT data not found: ' dftDataFile newline ...
        'Run kP_grid_analysis.m first.']);
end

load(dftDataFile);

% Frequencies up to this multiple of the excitation frequency are shown.
maximumDisplayedHarmonic = 5;

%% 2. SHOW THE FLOQUET-STABLE GRID

gridFigure = figure('Color','w');
gridAxes = axes(gridFigure);

stableAmplitudeIndices = find(any(stableOnSimulationGrid,2));
stableFrequencyIndices = find(any(stableOnSimulationGrid,1));

if isempty(stableAmplitudeIndices) || isempty(stableFrequencyIndices)
    error('No Floquet-stable points exist in the stored grid.');
end

displayedAmplitudeValues = amplitudeValues(stableAmplitudeIndices);
displayedFrequencyValues = frequencyValues(stableFrequencyIndices);
displayedStableRegion = ...
    stableOnSimulationGrid(stableAmplitudeIndices,stableFrequencyIndices);

imagesc(gridAxes, ...
    displayedFrequencyValues,displayedAmplitudeValues,displayedStableRegion);
axis(gridAxes,'xy');
xlabel(gridAxes,'excitation frequency, f [Hz]');
ylabel(gridAxes,'base displacement amplitude, A [m]');
title(gridAxes,{'Select a Floquet-stable point', ...
    'Click a green cell; press Enter to finish'});
colormap(gridAxes,[0.85 0.85 0.85; 0.20 0.65 0.30]);
caxis(gridAxes,[0 1]);
hold(gridAxes,'on');

selectionMarker = plot(gridAxes,NaN,NaN,'rx', ...
    'LineWidth',2,'MarkerSize',12);

%% 3. PREPARE THE DFT FIGURE

dftFigure = figure('Color','w');

amplitudeAxes = subplot(2,1,1);
amplitudeLine = plot(amplitudeAxes,NaN,NaN,'LineWidth',1.3);
grid(amplitudeAxes,'on');
xlabel(amplitudeAxes,'frequency [Hz]');
ylabel(amplitudeAxes,'angle amplitude [rad]');

phaseAxes = subplot(2,1,2);
phaseLine = plot(phaseAxes,NaN,NaN,'.-', ...
    'LineWidth',1.0,'MarkerSize',10);
grid(phaseAxes,'on');
xlabel(phaseAxes,'frequency [Hz]');
ylabel(phaseAxes,'phase [deg]');
ylim(phaseAxes,[-180 180]);

%% 4. SELECT AND DISPLAY POINTS

while ishandle(gridFigure)

    figure(gridFigure);
    [selectedFrequency,selectedAmplitude,mouseButton] = ginput(1);

    if isempty(mouseButton)
        break;
    end

    [~,amplitudeIndex] = min(abs(amplitudeValues-selectedAmplitude));
    [~,frequencyIndex] = min(abs(frequencyValues-selectedFrequency));

    actualAmplitude = amplitudeValues(amplitudeIndex);
    actualFrequency = frequencyValues(frequencyIndex);

    if ~stableOnSimulationGrid(amplitudeIndex,frequencyIndex)
        fprintf(['A = %.4f m, f = %.4f Hz is outside the ' ...
            'Floquet-stable region.\n'],actualAmplitude,actualFrequency);
        continue;
    end

    frequencyAxis = dftFrequencyHz(:,frequencyIndex);
    amplitudeSpectrum = dftAmplitude(:,amplitudeIndex,frequencyIndex);
    complexSpectrum = dftComplex(:,amplitudeIndex,frequencyIndex);
    phaseSpectrumDegrees = rad2deg(angle(complexSpectrum));

    % Phase is meaningless when its spectral amplitude is practically zero.
    phaseAmplitudeLimit = 1e-4*max(amplitudeSpectrum);
    phaseSpectrumDegrees(amplitudeSpectrum < phaseAmplitudeLimit) = NaN;

    maximumDisplayedFrequency = maximumDisplayedHarmonic*actualFrequency;
    displayedBins = frequencyAxis <= maximumDisplayedFrequency;

    set(selectionMarker,'XData',actualFrequency,'YData',actualAmplitude);

    set(amplitudeLine, ...
        'XData',frequencyAxis(displayedBins), ...
        'YData',amplitudeSpectrum(displayedBins));
    amplitudeTitle = sprintf( ...
        'DFT amplitude: A = %.4f m, f = %.4f Hz, theta0 = %.0f deg', ...
        actualAmplitude,actualFrequency,dftInitialAngleDegrees);
    title(amplitudeAxes,amplitudeTitle);
    xlim(amplitudeAxes,[0 maximumDisplayedFrequency]);

    set(phaseLine, ...
        'XData',frequencyAxis(displayedBins), ...
        'YData',phaseSpectrumDegrees(displayedBins));
    title(phaseAxes,'DFT phase relative to the excitation cosine');
    xlim(phaseAxes,[0 maximumDisplayedFrequency]);

    drawnow;
end
