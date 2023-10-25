function bc_plotLowAmpSpikes(chunkCentres, chunkLimits, spikeTimes, ...
    amplitudes, cutSpikes, means, STDs, thresholds, valid, thisUnit)

nBins = 20;
[~, binEdges] = histcounts(amplitudes, nBins);
binCentres = binEdges(1:end-1) + diff(binEdges(1:2))/2;
hists = NaN(size(chunkLimits,1), length(binCentres));
fits = NaN(size(chunkLimits,1), length(binCentres));

for ch = 1:size(chunkLimits,1)
    ind = spikeTimes >=chunkLimits(ch,1) & spikeTimes < chunkLimits(ch,2);
    amps = amplitudes(ind);
    amps(isnan(amps)) = [];
    hists(ch,:) = histcounts(amps, binEdges, 'Normalization', 'pdf');
    fits(ch,:) = bc_normpdf_cut(binCentres, means(ch), ...
        STDs(ch), thresholds(ch));
end

% if fitted means not within binEdges, pad histograms and fits with
% zeros to include fitted means
binDist = diff(binCentres(1:2));
if min(means) < binEdges(1)
    toAdd = ceil((binEdges(1) - min(means)) / binDist);
    binCentres = [binCentres(1)-toAdd*binDist : binDist: binCentres(1), ...
        binCentres(2:end)];
    binEdges = [binCentres, binCentres(end)+binDist] - binDist/2;
    hists = padarray(hists, [0 toAdd], 0, 'pre');
    fits = padarray(fits, [0 toAdd], 0, 'pre');
end
fits(~valid,:) = NaN;

mini = min([hists(:); fits(:)]);
maxi = max([hists(:); fits(:)]);

figure('WindowState', 'fullscreen')
tiledlayout(2,1)

nexttile
hold on
imagesc(mean(chunkCentres([1 end],:),2), binCentres([1 end]), hists', [mini maxi])
h = scatter(spikeTimes, amplitudes, 4, 'white', 'filled');
alpha(h, 0.3)
set(gca, 'YDir', 'normal')
xlim(chunkCentres([1 end]))
ylim(binEdges([1 end]))
colorbar
title('Histograms')

nexttile
hold on
imagesc(mean(chunkCentres([1 end],:),2), binCentres([1 end]), fits', [mini maxi])
h1 = errorbar(mean(chunkCentres,2), means, STDs, 'ro-', ...
    'MarkerFaceColor', 'r', 'LineWidth', 2, 'CapSize', 0);
h2 = scatter(spikeTimes(~cutSpikes), amplitudes(~cutSpikes), 4, 'white', 'filled');
alpha(h2, 0.3)
xlim(chunkCentres([1 end]))
ylim(binEdges([1 end]))
set(gca, 'YDir', 'normal')
colorbar
title('Fitted Gaussians')
xlabel('Time (s)')
ylabel('Amplitude')
legend(h1, 'Fitted means')

sgtitle(sprintf('Unit %d', thisUnit))