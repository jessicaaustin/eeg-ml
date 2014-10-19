function escaped = escape(S)

if iscell(S),
	escaped = cell(size(S));
	for i = 1:length(S),
		escaped{i} = regexprep(S{i}, '_', ' ');
	end
else
	escaped = regexprep(S, '_', ' ');
end
