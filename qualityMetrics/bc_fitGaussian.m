function [means, STDs, thresholds] = bc_fitGaussian(chunkLimits, spikeTimes, ...
    amplitudes)

means = NaN(size(chunkLimits,1),1);
STDs = NaN(size(chunkLimits,1),1);
thresholds = NaN(size(chunkLimits,1),1);
for ch = 1:size(chunkLimits,1)
    ind = spikeTimes >=chunkLimits(ch,1) & spikeTimes < chunkLimits(ch,2);
    amps = amplitudes(ind);
    amps(isnan(amps)) = [];
    thresholds(ch) = min(amps);
    gauss_pars = mle(amps, 'pdf', ...
        @(amps, mu, sigma) bc_normpdf_cut(amps, mu, sigma, thresholds(ch)), ...
        'Start', [mean(amps), std(amps)], ...
        'LowerBound', [0, 0], 'UpperBound', [max(amps), max(amps)]);
    means(ch) = gauss_pars(1);
    STDs(ch) = gauss_pars(2);
end