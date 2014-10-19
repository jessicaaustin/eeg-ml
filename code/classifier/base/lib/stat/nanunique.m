function U = nanunique(X)

I = ~isnan(X);
U = unique(X(I));
