function [chunkLimits, chunkCentres, invalidChunks, ampFilt, ampRanges, ...
    ampMins, ampThr] = ...
    bc_getOverlappingTimeChunks(start, stop, spikeTimes, amplitudes, param)

% TODO: introduce non-overlapping periods if spike amplitudes vary too fast
% 1. get smoothed time-dependent spike amplitude measure
% 2. take 1st derivative
% 3. introduce period borders if derivative crosses threshold
% 4. get overlapping chunks within each period

% Alternative
% 1. Estimate STD from amplitudes minus median filtered amplitudes
%    (possibly variable across time, but take care if amplitudes below detection threshold)
% 2. Divide these corrected amplitudes by STD

% 2. Subtract (param.maxPercSpikeMissing) * STD from median-filtered
% amplitudes -> time below 0 are invalid

% 2. Use median filtered amplitudes to detect times when STD is cut-off too
% much -> disregard
% 3. Use valid median filtered amplitudes to detect minimum -> use this
% point to cut off small spike amplitudes everywhere

% Or:
% Iterative process:
% caculate chunks and means as implemented
% subtract means and check how much median filtered amplitudes change
% if change is larger than 2 STD (median corrected amplitudes), then divide
% those chunks into smaller chunks


% parameters
centreSize = param.chunkCentreSize;
minSpikes = param.minNumSpikesPerChunk;
minSize = param.minChunkSize;
maxSize = param.maxChunkSize;
thresh = param.maxAmpChange;

% check that there are more spikes than minSpikes
if length(spikeTimes) < minSpikes
    chunkCentres = [start stop];
    chunkLimits = [start stop];
    invalidChunks = true;
    return
end

% estimate STD of straightened data
timeStep = 0.2; %in sec
spikeRate = length(spikeTimes) / (stop - start);
ampFilt = medfilt1(amplitudes, max(21, round(spikeRate * timeStep) * 4 + 1));
ampStraight = amplitudes - ampFilt;
ampSTD = std(ampStraight);
ampThr = ampSTD * thresh;

% define non-overlapping chunk centres
centreSize = (stop - start) / round((stop - start) / centreSize);
chunkCentres = (start : centreSize : stop)';
chunkCentres = [chunkCentres(1:end-1) chunkCentres(2:end)];
centreTimes = chunkCentres(:,1) + centreSize/2;

% find chunk for each non-overlapping centre with enough spikes and within
% min and max length
% 1. start by setting chunk size to min size
chunkLimits = centreTimes + [-minSize minSize]./2;
% 2. find number of spikes per chunk
numSpikes = NaN(size(chunkLimits,1), 1);
for ch = 1:size(chunkLimits,1)
    numSpikes(ch) = sum(spikeTimes >= chunkLimits(ch,1) & spikeTimes < chunkLimits(ch,2));
end
% 3. extend chunks to have at least minSpikes
tooSmall = find(numSpikes < minSpikes);
for ch = 1:length(tooSmall)
    sp = find(spikeTimes > centreTimes(tooSmall(ch)), 1);
    if isempty(sp)
        sp = length(spikeTimes);
    end
    sp_start = sp - ceil(minSpikes/2);
    sp_stop = sp_start + minSpikes - 1;
    if sp_start < 1
        sp_stop = sp_stop - sp_start + 1;
        sp_start = 1;
    end
    if sp_stop > length(spikeTimes)
        sp_start = sp_start - (sp_stop - length(spikeTimes));
        sp_stop = length(spikeTimes);
    end
    chunkLimits(tooSmall(ch),:) = [spikeTimes(sp_start) spikeTimes(sp_stop)+0.01];
end
% 4. check whether new sizes of chunks are larger maxSize or outside of
% chunkCentres
sizes = diff(chunkLimits, 1, 2);
invalidChunks = find(sizes > maxSize);
invalidChunks = unique([invalidChunks; find(chunkLimits(:,1) > chunkCentres(:,2) | ...
    chunkLimits(:,2) < chunkCentres(:,1))]);
% 5. check whether amplitudes in chunks show too much variability
ampRanges = NaN(size(chunkLimits,1), 1);
ampMins = NaN(size(chunkLimits,1), 1);
for ch = 1:size(chunkLimits,1)
    if ismember(ch, invalidChunks)
        continue
    end
    ind = spikeTimes >= chunkLimits(ch,1) & spikeTimes < chunkLimits(ch,2);
    ampRanges(ch) = range(ampFilt(ind));
    ampMins(ch) = min(ampFilt(ind));
end
invalidChunks = unique([invalidChunks; find(ampRanges > ampThr)]);