function [S] = chi_squared_sig_test(O, E, alpha);

S.values = O;
S.N = sum(S.values);
S.alpha = alpha;
S.test = 'chi_squared';

observed = O;
expected = E;
chi2stat = sum((observed-expected).^2 ./ expected);
S.P = 1 - chi2cdf(chi2stat, 1);

S.H = 0;
if (S.alpha > S.P),
    S.H = 1;
end