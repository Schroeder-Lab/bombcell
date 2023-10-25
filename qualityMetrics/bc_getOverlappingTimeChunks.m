function [chunkLimits, chunkCentres, invalidChunks] = ...
    bc_getOverlappingTimeChunks(start, stop, spikeTimes, param)

% parameters
centreSize = param.chunkCentreSize;
minSpikes = param.minNumSpikesPerChunk;
minSize = param.minChunkSize;
maxSize = param.maxChunkSize;

% check that there are more spikes than minSpikes
if length(spikeTimes) < minSpikes
    chunkCentres = [start stop];
    chunkLimits = [start stop];
    invalidChunks = true;
    return
end

% define non-overlapping chunk centres
centreSize = (stop - start) / round((stop - start) / centreSize);
chunkCentres = (start : centreSize : stop)';
chunkCentres = [chunkCentres(1:end-1) chunkCentres(2:end)];
centreTimes = chunkCentres(:,1) + centreSize/2;

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