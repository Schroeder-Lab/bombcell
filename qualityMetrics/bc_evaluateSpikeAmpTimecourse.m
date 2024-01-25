function [badAmpPeriods, lowAmpSpikes, goodAmpUnits] = ...
    bc_evaluateSpikeAmpTimecourse(param, ...
    spikeTimes_samples, spikeTemplates, templateAmplitudes, unitType, ...
    savePath, doPlot)
% SS
% ------
% Inputs
% ------
% param: parameter structure with fields:
%   tauR = 0.0010; %refractory period time (s)
%   tauC = 0.0002; %censored period time (s)
%   maxPercSpikesMissing: maximum percent (eg 30) of estimated spikes below detection
%       threshold to define timechunks in the recording on which to compute
%       quality metrics for each unit.
%   minNumSpikes: minimum number of spikes (eg 300) for unit to classify it as good
%   maxNtroughsPeaks: maximum number of troughs and peaks (eg 3) to classify unit
%       waveform as good
%   isSomatic: boolean, whether to keep only somatic spikes
%   maxRPVviolations: maximum estimated % (eg 20) of refractory period violations to classify unit as good
%   minAmplitude: minimum amplitude of raw waveform in microVolts to
%       classify unit as good
%   plotThis: boolean, whether to plot figures for each metric and unit - !
%       this will create * a lot * of plots if run on all units - use just
%       for debugging a particular issue / creating plots for one single
%       unit
%   rawFolder: string containing the location of the raw .dat or .bin file
%   deltaTimeChunk: size of time chunks to cut the recording in, in seconds
%       (eg 600 for 10 min time chunks or duration of recording if you don't
%       want time chunks)
%   ephys_sample_rate: recording sample rate (eg 30000)
%   nChannels: number of recorded channels, including any sync channels (eg
%       385)
%   nRawSpikesToExtract: number of spikes to extract from the raw data for
%       each waveform (eg 100)
%   nChannelsIsoDist: number of channels on which to compute the distance
%       metrics (eg 4)
%   computeDistanceMetrics: boolean, whether to compute distance metrics or not
%   isoDmin: minimum isolation distance to classify unit as single-unit
%   lratioMin: minimum l-ratio to classify unit as single-unit
%   ssMin: silhouette score to classify unit as single-unit
%   computeTimeChunks
%
% spikeTimes_samples: nSpikes × 1 uint64 vector giving each spike time in samples (*not* seconds)
%
% spikeTemplates: nSpikes × 1 uint32 vector giving the identity of each
%   spike's matched template
%
% templateAmplitudes: nSpikes × 1 double vector of the amplitude scaling factor
%   that was applied to the template when extracting that spike
%
% unitType: nUnits x 1 vector indicating whether each unit met the
%   threshold criterion to be classified as a single unit (1), noise
%   (0) or multi-unit (2)
%
%------
% Outputs
% ------
% badAmpPeriods: nPeriods x 3 double matrix; 1st column specifies unit
%   [1-nUnits]; 2nd and 3rd columns specify start and end (in seconds) of
%   "bad period" of that unit, i.e., period with too high spike amplitude
%   variation or period with too large cut off of low amplitude spikes.
%
% lowAmpSpikes: nSpikes x 1 logical vector classifying each spike as below
%   amplitude threshold (true) and above (false); spikes below threshold
%   should not be used for further analysis.
%
% goodAmpUnits: nUnits x 1 logical vector classifying all single units into
%   those with sufficient high amplitude spikes (true) or not enough high
%   amplitude spikes (false). Units not previously classified as single
%   unit (as in unitType) are NaN.

% get spike times
spikeTimes_seconds = spikeTimes_samples ./ param.ephys_sample_rate;

% get unique templates 
uniqueTemplates = unique(spikeTemplates);

badAmpPeriods = [];
lowAmpSpikes = false(size(spikeTimes_seconds));
goodAmpUnits = NaN(length(unitType), 1);

for iUnit = 1:length(unitType)
    if unitType(iUnit) ~= 1
        continue
    end
    clearvars thisUnit theseSpikeTimes theseAmplis theseSpikeTemplates

    % get this unit's attributes 
    thisUnit = uniqueTemplates(iUnit);
    theseSpikeInds = find(spikeTemplates == thisUnit);
    theseSpikeTimes = spikeTimes_seconds(theseSpikeInds);
    theseAmplis = templateAmplitudes(theseSpikeInds);
    % sort spike times in ascending order
    [theseSpikeTimes, order] = sort(theseSpikeTimes);
    theseSpikeInds = theseSpikeInds(order);
    theseAmplis = theseAmplis(order);

    %% cut off low amplitude spikes to create consistent percentage spikes missing across time
    % determine overlapping time chunks
    % Note: invalidChunks are those with too few spikes, or with too much
    % variation in spike amplitudes
    [chunkLimits, chunkCentres, invalidChunks, ampFilt, ampRanges, ampMins, ampThr] = ...
        bc_getOverlappingTimeChunks(min(spikeTimes_seconds), ...
        max(spikeTimes_seconds), theseSpikeTimes, theseAmplis, param);

    % for each time chunk: fit Gaussian, get cut-off in terms of STDs
    valid = true(size(chunkLimits,1), 1);
    valid(invalidChunks) = false;
    chunkAmpMeans = NaN(size(valid));
    chunkAmpSTDs = NaN(size(valid));
    chunkAmpCutOffs = NaN(size(valid));
    [chunkAmpMeans(valid), chunkAmpSTDs(valid), chunkAmpCutOffs(valid)] = ...
        bc_fitGaussian(chunkLimits(valid,:), theseSpikeTimes, theseAmplis);

    % get percentage missing spikes
    cutOffs = normcdf(chunkAmpCutOffs, chunkAmpMeans, chunkAmpSTDs) .* 100;

    % disregard chunks with missing spikes > maxPercSpikesMissing
    valid(cutOffs > param.maxPercSpikesMissing) = false;

    % determine whether unit is valid
    if sum(valid) / length(valid) >= param.minValid
        goodAmpUnits(iUnit) = 1;
    else
        goodAmpUnits(iUnit) = 0;
    end

    % concatenate consecutive invalid chunks, return start and end of each
    onoff = diff(valid);
    if onoff(find(onoff,1)) == 1
        onoff(1) = -1;
    end
    if onoff(find(onoff, 1, 'last')) == -1
        onoff(end) = 1;
    end
    onsets = find(onoff == -1);
    offsets = find(onoff == 1);
    for ind = 1:length(onsets)
        badAmpPeriods(end+1,:) = [thisUnit, ...
            chunkCentres(onsets(ind),1), chunkCentres(offsets(ind),2)];
    end

    % find low amplitude spikes that should be cut off to guarantee the
    % same cut-off across time
    if all(~valid)
        cutSpikes = true(size(theseSpikeTimes));
    else
        cutSpikes = false(size(theseSpikeTimes));

        % determine max cut off
        maxCutOff = max(cutOffs(valid));
        maxCutOff = norminv(maxCutOff/100); % in terms of STDs from mean

        % for each time chunk: retrieve spike IDs with amplitudes smaller than
        % max STD cut-off
        cutOffs_amps = chunkAmpMeans + maxCutOff .* chunkAmpSTDs;
        for ch = 1:length(chunkCentres)
            ind = theseSpikeTimes >= chunkCentres(ch,1) & ...
                theseSpikeTimes < chunkCentres(ch,2);
            if valid(ch)
                cutSpikes = cutSpikes | (ind & theseAmplis < cutOffs_amps(ch));
            else
                cutSpikes = cutSpikes | ind;
            end
        end
    end
    lowAmpSpikes(theseSpikeInds(cutSpikes)) = true;

    if doPlot
        % plot old and new data (including Gaussians)
        bc_plotLowAmpSpikes(chunkCentres, chunkLimits, theseSpikeTimes, ...
            theseAmplis, ampFilt, ampRanges, ampMins, ampThr, cutSpikes, ...
            chunkAmpMeans, chunkAmpSTDs, chunkAmpCutOffs, valid, thisUnit);
    end
end

% save data
writeNPY(badAmpPeriods, fullfile(savePath, 'badAmplitudeIntervals.npy'));
writeNPY(lowAmpSpikes, fullfile(savePth, 'lowAmplitudeSpikes.npy'));
writeNPY(goodAmpUnits, fullfile(savePath, 'goodAmplitudeUnits.npy'));