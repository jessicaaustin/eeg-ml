function [f, u] = nanfreq(x)

I = ~isnan(x);
[f u] = freq(x(I));
