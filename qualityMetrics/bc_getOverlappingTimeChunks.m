function [chunkLimits, centreLimits, invalidChunks] = ...
    bc_getOverlappingTimeChunks(start, stop, spikeTimes, params)

% parameters
centreSize = params.chunkCentreSize;
minSpikes = params.minNumSpikesPerChunk;
minSize = params.minChunkSize;
maxSize = params.maxChunkSize;

% check that there are more spikes than minSpikes
if length(spikeTimes) < minSpikes
    centreLimits = NaN;
    chunkLimits = NaN;
    return
end

% define non-overlapping chunk centres
centreLimits = (start : centreSize : stop)';
centreLimits = [centreLimits(1:end-1) centreLimits(2:end)];
centreTimes = centreLimits(:,1) + centreSize/2;

% find chunk for each non-overlapping centre with enough spikes and within
% min and max length
% 1. start by setting chunk size to min size
chunkLimits = centreTimes + [-minSize minSize]./2;
% 2. find chunks with fewer spikes than minSpikes
numSpikes = NaN(size(chunkLimits,1), 1);
for ch = 1:size(chunkLimits,1)
    numSpikes(ch) = sum(spikeTimes >=chunkLimits(ch,1) & spikeTimes < chunkLimits(ch,2));
end
% 3. extend chunks to have at least minSpikes
tooSmall = find(numSpikes < minSpikes);
for ch = 1:length(tooSmall)
    sp = find(spikeTimes > centreTimes(tooSmall(ch)), 1);
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
% 4. check whether new sizes of chunks are larger maxSize
sizes = diff(chunkLimits, 1, 2);
invalidChunks = find(sizes > maxSize);