function LIA = issubstr(A, B)

LIA = [];
for b = 1:length(B),
	I = strfind(A, B{b});
	for i = 1:length(I),
		if ~isempty(I{i}),
			LIA(end+1) = i;
		end
	end % i
end % b
