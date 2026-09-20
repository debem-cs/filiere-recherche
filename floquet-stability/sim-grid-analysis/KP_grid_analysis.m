%% KAPITZA PENDULUM: ANALYSE STORED GRID RESULTS
%
% This script does not run KP_sim.slx. It only loads the existing MAT files,
% calculates derived quantities and creates figures in the generated folder.

clear;
clc;
close all;

%% 1. LOAD THE EXISTING RESULTS

projectFolder = fileparts(mfilename('fullpath'));
outputFolder = fullfile(projectFolder,'generated');

gridResultsFile = fullfile(outputFolder,'KP_grid_results.mat');
dftResultsFile = fullfile(outputFolder,'KP_DFT_results.mat');
stabilityResultsFile = fullfile(outputFolder,'kapitza_stability_data.mat');

if ~isfile(gridResultsFile)
    error(['Grid results not found: ' gridResultsFile]);
end

if ~isfile(dftResultsFile)
    error(['DFT results not found: ' dftResultsFile]);
end

if ~isfile(stabilityResultsFile)
    error(['Stability results not found: ' stabilityResultsFile]);
end

gridData = load(gridResultsFile);
dftData = load(dftResultsFile);
stabilityData = load(stabilityResultsFile);

amplitudeValues = gridData.amplitudeValues;
frequencyValues = gridData.frequencyValues;
initialAngleDegrees = gridData.initialAngleDegrees;
referenceSettlingAngleDegrees = gridData.referenceSettlingAngleDegrees;
settlingTime = gridData.settlingTime;
settlingTimeAllAngles = gridData.settlingTimeAllAngles;
converged = gridData.converged;
attractionAngle = gridData.attractionAngle;
stableOnSimulationGrid = gridData.stableOnSimulationGrid;

dftAmplitude = dftData.dftAmplitude;
dftComplex = dftData.dftComplex;
dftFrequencyHz = dftData.dftFrequencyHz;

if isfield(dftData,'dftInitialAngleDegrees')
    dftInitialAngleDegrees = dftData.dftInitialAngleDegrees;
else
    dftInitialAngleDegrees = 5;
end

if ~isequal(amplitudeValues,dftData.amplitudeValues) || ...
        ~isequal(frequencyValues,dftData.frequencyValues)
    error('KP_grid_results.mat and KP_DFT_results.mat use different grids.');
end

%% 2. CROP THE FIGURES TO THE STABLE INTERVAL

stableAmplitudeIndices = find(any(stableOnSimulationGrid,2));
stableFrequencyIndices = find(any(stableOnSimulationGrid,1));

if isempty(stableAmplitudeIndices) || isempty(stableFrequencyIndices)
    error('No Floquet-stable points exist in the stored grid.');
end


displayedAmplitudeValues = amplitudeValues(stableAmplitudeIndices);
displayedFrequencyValues = frequencyValues(stableFrequencyIndices);
displayedStableRegion = ...
    stableOnSimulationGrid(stableAmplitudeIndices,stableFrequencyIndices);

%% 3. ORIGINAL SETTLING-TIME FIGURE

displayedSettlingTime = ...
    settlingTime(stableAmplitudeIndices,stableFrequencyIndices);

settlingFigure = figure('Color','w');
settlingImage = imagesc( ...
    displayedFrequencyValues,displayedAmplitudeValues,displayedSettlingTime);
set(settlingImage,'AlphaData',~isnan(displayedSettlingTime));
axis xy;
set(gca,'Color',[0.85 0.85 0.85]);
xlabel('excitation frequency, f [Hz]');
ylabel('base displacement amplitude, A [m]');
title(sprintf('Settling time from \\theta(0) = %.0f deg [s]', ...
    referenceSettlingAngleDegrees));
colorbar;

exportgraphics(settlingFigure, ...
    fullfile(outputFolder,'KP_settling_time.png'),'Resolution',200);

%% 4. ORIGINAL ATTRACTION-ANGLE FIGURE

displayedAttractionAngle = ...
    attractionAngle(stableAmplitudeIndices,stableFrequencyIndices);

attractionFigure = figure('Color','w');
attractionImage = imagesc( ...
    displayedFrequencyValues,displayedAmplitudeValues,displayedAttractionAngle);
set(attractionImage,'AlphaData',displayedStableRegion);
axis xy;
set(gca,'Color',[0.85 0.85 0.85]);
xlabel('excitation frequency, f [Hz]');
ylabel('base displacement amplitude, A [m]');
title('Largest tested initial angle that converges [deg]');
colorbar;
caxis([0 max(initialAngleDegrees)]);

exportgraphics(attractionFigure, ...
    fullfile(outputFolder,'KP_attraction_region.png'),'Resolution',200);

%% 5. ANALYSIS SETTINGS

% The strongest spectral peak a1 is always relevant. Every additional peak
% must be greater than 10% of a1.
relativePeakThreshold = 0.10;
minimumPeakSeparationBins = 2;
maximumAnalyzedHarmonic = 5;
numberOfPeaksStored = 3;

numberOfAmplitudes = length(amplitudeValues);
numberOfFrequencies = length(frequencyValues);

systemOscillationEnergyPerMass = ...
    NaN(numberOfAmplitudes,numberOfFrequencies);
controlSignalEnergyPerPeriod = ...
    NaN(numberOfAmplitudes,numberOfFrequencies);
relevantFrequencyCount = ...
    NaN(numberOfAmplitudes,numberOfFrequencies);

relevantFrequencyHz = ...
    NaN(numberOfAmplitudes,numberOfFrequencies,numberOfPeaksStored);
relevantFrequencyRatio = ...
    NaN(numberOfAmplitudes,numberOfFrequencies,numberOfPeaksStored);
relevantPeakAmplitude = ...
    NaN(numberOfAmplitudes,numberOfFrequencies,numberOfPeaksStored);
relativePeakAmplitude = ...
    NaN(numberOfAmplitudes,numberOfFrequencies,numberOfPeaksStored);

%% 6. ENERGY AND RELEVANT-FREQUENCY CALCULATIONS

for amplitudeIndex = 1:numberOfAmplitudes
    for frequencyIndex = 1:numberOfFrequencies

        if ~stableOnSimulationGrid(amplitudeIndex,frequencyIndex)
            continue;
        end

        excitationFrequency = frequencyValues(frequencyIndex);
        excitationAngularFrequency = 2*pi*excitationFrequency;
        baseAmplitude = amplitudeValues(amplitudeIndex);

        frequencyAxis = dftFrequencyHz(:,frequencyIndex);
        amplitudeSpectrum = dftAmplitude(:,amplitudeIndex,frequencyIndex);

        if any(isnan(amplitudeSpectrum)) || any(isnan(frequencyAxis))
            continue;
        end

        % Parseval calculation using the one-sided amplitude spectrum.
        % Nonzero, non-Nyquist sinusoidal components contribute A_k^2/2.
        parsevalWeight = 0.5*ones(size(amplitudeSpectrum));
        parsevalWeight(1) = 1;
        parsevalWeight(end) = 1;

        meanSquareAngle = sum( ...
            parsevalWeight.*amplitudeSpectrum.^2);
        meanSquareAngularVelocity = sum( ...
            parsevalWeight.*(2*pi*frequencyAxis.*amplitudeSpectrum).^2);

        systemOscillationEnergyPerMass(amplitudeIndex,frequencyIndex) = ...
            0.5*stabilityData.L^2*meanSquareAngularVelocity ...
            + 0.5*stabilityData.g*stabilityData.L*meanSquareAngle;

        % Signal energy of u(t)=a0*cos(omega*t) during one forcing period.
        accelerationAmplitude = ...
            baseAmplitude*excitationAngularFrequency^2;
        excitationPeriod = 1/excitationFrequency;
        controlSignalEnergyPerPeriod(amplitudeIndex,frequencyIndex) = ...
            0.5*accelerationAmplitude^2*excitationPeriod;

        % Ignore DC and frequencies above the selected harmonic limit.
        analyzedBins = find( ...
            frequencyAxis > 0 & ...
            frequencyAxis <= maximumAnalyzedHarmonic*excitationFrequency);

        analyzedAmplitude = amplitudeSpectrum(analyzedBins);

        if isempty(analyzedAmplitude)
            relevantFrequencyCount(amplitudeIndex,frequencyIndex) = 0;
            continue;
        end

        candidateBins = [];
        candidateAmplitudes = [];

        for analyzedIndex = 1:length(analyzedBins)

            spectrumIndex = analyzedBins(analyzedIndex);

            if spectrumIndex == 1 || spectrumIndex == length(amplitudeSpectrum)
                continue;
            end

            presentAmplitude = amplitudeSpectrum(spectrumIndex);

            isLocalMaximum = ...
                presentAmplitude > amplitudeSpectrum(spectrumIndex-1) && ...
                presentAmplitude >= amplitudeSpectrum(spectrumIndex+1);

            if isLocalMaximum
                candidateBins(end+1) = spectrumIndex; %#ok<SAGROW>
                candidateAmplitudes(end+1) = presentAmplitude; %#ok<SAGROW>
            end
        end

        % Always include the largest spectral bin. This also guarantees a1
        % when the largest value occurs at the edge of the analyzed range.
        [largestAmplitude,largestAmplitudeIndex] = max(analyzedAmplitude);
        largestAmplitudeBin = analyzedBins(largestAmplitudeIndex);

        if ~any(candidateBins == largestAmplitudeBin)
            candidateBins(end+1) = largestAmplitudeBin;
            candidateAmplitudes(end+1) = largestAmplitude;
        end

        [candidateAmplitudes,candidateOrder] = ...
            sort(candidateAmplitudes,'descend');
        candidateBins = candidateBins(candidateOrder);

        % a1 is always retained. All other candidates are compared with a1.
        relevanceLimit = relativePeakThreshold*candidateAmplitudes(1);
        relevantCandidates = false(size(candidateAmplitudes));
        relevantCandidates(1) = true;

        if length(candidateAmplitudes) > 1
            relevantCandidates(2:end) = ...
                candidateAmplitudes(2:end) > relevanceLimit;
        end

        candidateBins = candidateBins(relevantCandidates);
        candidateAmplitudes = candidateAmplitudes(relevantCandidates);

        selectedBins = [];
        selectedAmplitudes = [];

        for candidateIndex = 1:length(candidateBins)

            presentBin = candidateBins(candidateIndex);

            if isempty(selectedBins) || all( ...
                    abs(presentBin-selectedBins) > minimumPeakSeparationBins)
                selectedBins(end+1) = presentBin; %#ok<SAGROW>
                selectedAmplitudes(end+1) = ...
                    candidateAmplitudes(candidateIndex); %#ok<SAGROW>
            end
        end

        relevantFrequencyCount(amplitudeIndex,frequencyIndex) = ...
            length(selectedBins);

        peaksToStore = min(numberOfPeaksStored,length(selectedBins));

        for peakIndex = 1:peaksToStore

            presentBin = selectedBins(peakIndex);
            presentFrequency = frequencyAxis(presentBin);
            presentAmplitude = selectedAmplitudes(peakIndex);

            relevantFrequencyHz(amplitudeIndex,frequencyIndex,peakIndex) = ...
                presentFrequency;
            relevantFrequencyRatio(amplitudeIndex,frequencyIndex,peakIndex) = ...
                presentFrequency/excitationFrequency;
            relevantPeakAmplitude(amplitudeIndex,frequencyIndex,peakIndex) = ...
                presentAmplitude;

            if peakIndex == 1
                relativePeakAmplitude(amplitudeIndex,frequencyIndex,peakIndex) = 1;
            elseif selectedAmplitudes(1) > 0
                relativePeakAmplitude(amplitudeIndex,frequencyIndex,peakIndex) = ...
                    presentAmplitude/selectedAmplitudes(1);
            end
        end
    end
end

%% 7. RESIDUAL SYSTEM-OSCILLATION ENERGY MAP

displayedSystemEnergy = systemOscillationEnergyPerMass( ...
    stableAmplitudeIndices,stableFrequencyIndices);
displayedLogSystemEnergy = NaN(size(displayedSystemEnergy));
positiveSystemEnergy = displayedSystemEnergy > 0;
displayedLogSystemEnergy(positiveSystemEnergy) = ...
    log10(displayedSystemEnergy(positiveSystemEnergy));

systemEnergyFigure = figure('Color','w');
systemEnergyImage = imagesc(displayedFrequencyValues, ...
    displayedAmplitudeValues,displayedLogSystemEnergy);
set(systemEnergyImage,'AlphaData', ...
    displayedStableRegion & ~isnan(displayedLogSystemEnergy));
axis xy;
set(gca,'Color',[0.85 0.85 0.85]);
xlabel('excitation frequency, f [Hz]');
ylabel('base displacement amplitude, A [m]');
title(sprintf( ...
    'Residual oscillation energy, final DFT from \\theta(0) = %.0f deg', ...
    dftInitialAngleDegrees));
systemEnergyColorbar = colorbar;
systemEnergyColorbar.Label.String = ...
    'log_{10} specific oscillation energy [J/kg]';

exportgraphics(systemEnergyFigure, ...
    fullfile(outputFolder,'KP_system_oscillation_energy.png'),'Resolution',200);

%% 8. CONTROL-SIGNAL ENERGY MAP

displayedControlEnergy = controlSignalEnergyPerPeriod( ...
    stableAmplitudeIndices,stableFrequencyIndices);
displayedLogControlEnergy = log10(displayedControlEnergy);

controlEnergyFigure = figure('Color','w');
controlEnergyImage = imagesc(displayedFrequencyValues, ...
    displayedAmplitudeValues,displayedLogControlEnergy);
set(controlEnergyImage,'AlphaData',displayedStableRegion);
axis xy;
set(gca,'Color',[0.85 0.85 0.85]);
xlabel('excitation frequency, f [Hz]');
ylabel('base displacement amplitude, A [m]');
title('Control acceleration-signal energy per period');
controlEnergyColorbar = colorbar;
controlEnergyColorbar.Label.String = ...
    'log_{10} integral of u^2 dt [(m/s^2)^2 s]';

exportgraphics(controlEnergyFigure, ...
    fullfile(outputFolder,'KP_control_signal_energy.png'),'Resolution',200);

%% 9. NUMBER OF RELEVANT FREQUENCIES

displayedFrequencyCount = relevantFrequencyCount( ...
    stableAmplitudeIndices,stableFrequencyIndices);

frequencyCountFigure = figure('Color','w');
frequencyCountImage = imagesc(displayedFrequencyValues, ...
    displayedAmplitudeValues,displayedFrequencyCount);
set(frequencyCountImage,'AlphaData', ...
    displayedStableRegion & ~isnan(displayedFrequencyCount));
axis xy;
set(gca,'Color',[0.85 0.85 0.85]);
xlabel('excitation frequency, f [Hz]');
ylabel('base displacement amplitude, A [m]');
title('Number of relevant oscillation frequencies');
frequencyCountColorbar = colorbar;
frequencyCountColorbar.Label.String = 'number of spectral peaks';
storedPeakCounts = relevantFrequencyCount(stableOnSimulationGrid);
storedPeakCounts = storedPeakCounts(~isnan(storedPeakCounts));

if isempty(storedPeakCounts)
    maximumPeakCount = 1;
else
    maximumPeakCount = max(storedPeakCounts);
end

caxis([0 max(1,maximumPeakCount)]);

exportgraphics(frequencyCountFigure, ...
    fullfile(outputFolder,'KP_relevant_frequency_count.png'),'Resolution',200);

%% 10. THREE STRONGEST RELEVANT FREQUENCIES

frequencyRatioFigure = figure('Color','w');
frequencyRatioFigure.Position(3:4) = [1450 480];

storedFrequencyRatios = relevantFrequencyRatio(~isnan(relevantFrequencyRatio));

if isempty(storedFrequencyRatios)
    frequencyRatioColorLimit = 1;
else
    frequencyRatioColorLimit = max(1,ceil(max(storedFrequencyRatios)));
end

for peakIndex = 1:numberOfPeaksStored

    presentFrequencyRatio = relevantFrequencyRatio( ...
        stableAmplitudeIndices,stableFrequencyIndices,peakIndex);

    presentAxes = subplot(1,numberOfPeaksStored,peakIndex);
    presentImage = imagesc(presentAxes,displayedFrequencyValues, ...
        displayedAmplitudeValues,presentFrequencyRatio);
    set(presentImage,'AlphaData',~isnan(presentFrequencyRatio));
    axis(presentAxes,'xy');
    set(presentAxes,'Color',[0.85 0.85 0.85]);
    xlabel(presentAxes,'excitation frequency, f [Hz]');
    ylabel(presentAxes,'base displacement amplitude, A [m]');
    title(presentAxes,sprintf('Peak %d: frequency ratio f_%d/f', ...
        peakIndex,peakIndex));
    presentColorbar = colorbar(presentAxes);
    presentColorbar.Label.String = sprintf('f_%d/f',peakIndex);
    caxis(presentAxes,[0 frequencyRatioColorLimit]);
end

colormap(frequencyRatioFigure,turbo(256));

exportgraphics(frequencyRatioFigure, ...
    fullfile(outputFolder,'KP_relevant_frequency_ratios.png'),'Resolution',200);

%% 11. STRENGTH OF THE SECOND AND THIRD PEAKS

secondaryStrengthFigure = figure('Color','w');
secondaryStrengthFigure.Position(3:4) = [1050 480];

for peakIndex = 2:numberOfPeaksStored

    presentRelativeAmplitude = relativePeakAmplitude( ...
        stableAmplitudeIndices,stableFrequencyIndices,peakIndex);

    presentAxes = subplot(1,numberOfPeaksStored-1,peakIndex-1);
    presentImage = imagesc(presentAxes,displayedFrequencyValues, ...
        displayedAmplitudeValues,presentRelativeAmplitude);
    set(presentImage,'AlphaData',~isnan(presentRelativeAmplitude));
    axis(presentAxes,'xy');
    set(presentAxes,'Color',[0.85 0.85 0.85]);
    xlabel(presentAxes,'excitation frequency, f [Hz]');
    ylabel(presentAxes,'base displacement amplitude, A [m]');
    title(presentAxes,sprintf('Relative strength of peak %d',peakIndex));
    presentColorbar = colorbar(presentAxes);
    presentColorbar.Label.String = sprintf('a_%d/a_1',peakIndex);
    caxis(presentAxes,[0 1]);
end

colormap(secondaryStrengthFigure,turbo(256));

exportgraphics(secondaryStrengthFigure, ...
    fullfile(outputFolder,'KP_secondary_peak_strength.png'),'Resolution',200);

%% 12. SAVE THE DERIVED ANALYSIS RESULTS

save(fullfile(outputFolder,'KP_analysis_results.mat'), ...
    'systemOscillationEnergyPerMass', ...
    'controlSignalEnergyPerPeriod', ...
    'relevantFrequencyCount', ...
    'relevantFrequencyHz', ...
    'relevantFrequencyRatio', ...
    'relevantPeakAmplitude', ...
    'relativePeakAmplitude', ...
    'relativePeakThreshold', ...
    'minimumPeakSeparationBins', ...
    'maximumAnalyzedHarmonic', ...
    'frequencyRatioColorLimit', ...
    'numberOfPeaksStored');
