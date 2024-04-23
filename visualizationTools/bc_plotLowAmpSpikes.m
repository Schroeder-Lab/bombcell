function bc_plotLowAmpSpikes(chunkCentres, chunkLimits, spikeTimes, ...
    amplitudes, ampFilt, ampRanges, ampMins, ampThr, cutSpikes, means, STDs, thresholds, ...
    valid, thisUnit)

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
t = mean(chunkCentres,2);

figure('WindowState', 'maximized')
tiledlayout(2,1)
ax = zeros(1,2);

ax(1) = nexttile;
hold on
imagesc(t, binCentres([1 end]), hists', [mini maxi])
h = scatter(spikeTimes, amplitudes, 4, 'white', 'filled');
alpha(h, 0.7)
set(gca, 'YDir', 'normal')
xlim(chunkCentres([1 end]))
ylim(binEdges([1 end]))
c = colorbar;
c.Label.String = 'Measured frequency';
title('Histograms')

ax(2) = nexttile;
hold on
imagesc(t, binCentres([1 end]), fits', [mini maxi])
h1 = errorbar(t, means, STDs, 'ro-', ...
    'MarkerFaceColor', 'r', 'LineWidth', 2, 'CapSize', 0);
h2 = plot(mean(chunkCentres(~valid,:),2), means(~valid), 'wx', 'LineWidth', 2);
h3 = scatter(spikeTimes(~cutSpikes), amplitudes(~cutSpikes), 4, 'white', 'filled');
h4 = plot(spikeTimes, ampFilt, 'w', 'LineWidth', 0.5);
h5 = plot(t, ampMins, 'w', 'LineWidth', 2);
     plot(t, ampMins + ampRanges, 'w', 'LineWidth', 2);
h6 = plot(t, ampMins + ampThr, 'w-.', 'LineWidth', 1);
alpha(h3, 0.3)
xlim(chunkCentres([1 end]))
ylim(binEdges([1 end]))
set(gca, 'YDir', 'normal')
c = colorbar;
c.Label.String = 'Fitted frequency';
title('Fitted Gaussians')
xlabel('Time (s)')
ylabel('Amplitude')
if all(valid)
    l = legend([h1 h4 h5 h6], 'Gaussian', 'Median', 'Range', 'Threshold');
else
    l = legend([h1 h2 h4 h5 h6], 'Gaussian', 'Cut off', 'Median', 'Range', 'Threshold');
end
l.Color = 'none';
l.TextColor = 'w';
l.Box = 'off';

sgtitle(sprintf('Unit %d', thisUnit))
linkaxes(ax, 'x')