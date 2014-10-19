function D = nanpdist(X, DISTANCE)

[m n] = size(X);

c = 0;
for i = 1:m,
	for j = i+1:m,
		c = c + 1;

		I = find(and(~isnan(X(i,:)), ~isnan(X(j,:))));

		if length(unique(X(i,I))) < 2 || length(unique(X(j,I))) < 2,
			D(c) = NaN;
		else,
			D(c) = pdist([X(i,I); X(j,I)], DISTANCE);
		end
	end % j
end % i
