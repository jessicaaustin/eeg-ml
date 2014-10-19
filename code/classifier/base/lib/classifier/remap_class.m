function Y = remap_class(X, newclass)

Y = NaN(length(X), 1);

for i = 1:length(X),
	for j = 1:length(newclass),
		if strlocate(newclass{j}, X{t}) ~= -1,
			Y(i) = j;
		end
	end % j
end % i
