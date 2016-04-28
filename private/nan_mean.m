function out = nan_mean(in)

nans = find(isnan(in));
in(nans) = 0;
sums = sum(in);
nonnans = ones(size(in));
nonnans(nans) = 0;
nonnans = sum(nonnans);
nononnans = find(nonnans==0);
nonnans(nononnans) = 1;
out = sum(in)./nonnans;
out(nononnans) = NaN;


