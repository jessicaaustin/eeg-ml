for h = keys(H),
	k = h{:};
	k = strrep(k, '+', '_');
	k = strrep(k, '-', '_');
	k = strrep(k, '*', '_');
	k = strrep(k, '/', '_');

	eval( sprintf('%s = H(''%s'');', k, h{:}) );
	%eval( sprintf('%s = H(''%s'');', h{:}, h{:}) );
end
