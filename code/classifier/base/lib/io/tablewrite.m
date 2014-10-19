function tablewrite(file, M, header, rows)

fid = fopen(file, 'w');

if nargin > 2,
	%fprintf(fid, '\t');
	for j = 1:size(M,2),
		fprintf(fid, '%s\t', header{j});
	end
	fprintf(fid, '\n');
end

for i = 1:size(M,1),
	if nargin > 3,
		fprintf(fid, '%s\t', rows{i});
	end
	for j = 1:size(M,2),
		fprintf(fid, '%6d\t', M(i,j));
	end
	fprintf(fid, '\n');
end

fclose(fid);
