function [f, u] = freq(x)

if isempty(x), f = []; u = []; return; end

x = reshape(x, 1, []);

u = unique(x);

if iscell(x) && ischar(x{1}),
	for i = 1:length(u),
		f(i) = sum(strcmp(x, u(i)));
	end
else,
	for i = 1:length(u),
		f(i) = sum(x == u(i));
	end
end
