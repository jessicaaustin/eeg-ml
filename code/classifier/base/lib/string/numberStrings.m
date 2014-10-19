function S = numberStrings(S)

for i = 1:length(S),
	S{i} = sprintf('%d %s', i, S{i});
end
