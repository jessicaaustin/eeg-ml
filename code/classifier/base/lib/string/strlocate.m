function I = strlocate(S, str)

results = find(strcmp(S, str), 1);

if isempty(results),
	I = -1;
else,
	I = results;
end
