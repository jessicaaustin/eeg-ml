for li = 1:length(list),
	k = list{li};
	k = strrep(k, '+', '_');
	k = strrep(k, '-', '_');
	k = strrep(k, '*', '_');
	k = strrep(k, '/', '_');
	eval(sprintf('%s = %d;', k, (enum - 1) + li));
end