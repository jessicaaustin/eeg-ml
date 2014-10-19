function varargout = sscanp(s, p)

if strcmp(substr(version, 0,1), '6'),
	[start finish tokens] = regexp(s, p);
	for i = 1:size(tokens{:},1),
		varargout{i} = s(tokens{:}(i,1):tokens{:}(i,2));
	end
else,
	[i1 i2 extents matches tokens] = regexp(s, p);
	for i = 1:length(tokens{:}),
		varargout{i} = tokens{:}{i};
	end
end
