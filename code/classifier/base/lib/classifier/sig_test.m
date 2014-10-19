function [S] = sig_test(X, chance, alpha, tail);

S.values = X;
S.mean = nanmean(X);
S.std = nanstd(X);
S.N = sum(~isnan(X));

if nargin > 1,

S.test = 'ttest';
S.chance = chance;
S.alpha = alpha;
S.tail = tail;

[S.H S.P] = ttest(X, chance, alpha, tail);
	if isnan(S.H), S.H = 0; end;

end
