function y = bc_normpdf_cut(x, mu, sigma, threshold)

y = NaN(size(x));
y(x < threshold) = 0;
y(x >= threshold) = normpdf(x(x >= threshold), mu, sigma);
cut_area = normcdf(threshold, mu, sigma);
y = y ./ (1 - cut_area);

end