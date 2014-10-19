function tablewrite2(filename, H, M, T);

fid = fopen(filename, 'w');

for h = 1:length(H),
	fprintf(fid, '%s\t', H{h});
end % h

fprintf(fid, '\n');

for t = 1:size(M,1),
	for h = 1:length(H),
		if isfield(T, H{h}),
			fprintf(fid, '%s\t', eval(sprintf('T.%s{t}', H{h})));
		else,
			fprintf(fid, '%.4f\t', M(t,h));
		end
	end % h
	fprintf(fid, '\n');
end % t

fclose(fid);
