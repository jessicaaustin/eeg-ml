function [I2] = balance_class(Y, I, mode)

I2 = [];

[f u] = freq(Y(I));

for i = 1:length(u),
	Idx = find(and(I, Y == u(i)));

	switch mode,
	case 'oversample',
		N = max(f);
		if f(i) < N,
			r = randsample(f(i), N, true);
		else,
			r = 1:f(i);
		end

	case 'undersample',
		N = min(f);
		if f(i) > N,
			r = randsample(f(i), N, false);
		else,
			r = 1:f(i);
		end

	case 'truncate',
		N = min(f);
		if f(i) > N,
			r = 1:N;
		else,
			r = 1:f(i);
		end

	otherwise,
        I2 = I;
        return;

	end % switch

	I2  = [I2; Idx(r)];
end % i
