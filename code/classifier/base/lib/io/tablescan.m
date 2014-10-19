function [M, H, C, R, S] = tablescan(file, varargin)

% M is table converted to a numeric matrix
% H is header of the table
% C is a column view of the table
% R is a row view of the table
% S is a structure view of the table

delimiter = '	';

args = varargin;
for i = 1:2:length(args),
	eval(sprintf('%s = args{i+1};', args{i}));
end

lc = linecount(file);

fid = fopen(file);

% 1. H

line = fgetl(fid);
headers = tokenizer(line, delimiter);

	% Some fix to enable reading NeuroView files
	for i = 1:length(headers),
		if strcmp(headers{i}(1), '%'),
			headers{i} = headers{i}(2:end);
		elseif ~isnan(str2double(headers{i})),
			headers{i} = sprintf('H%s', strrep(headers{i}, '.', '_'));
		end
	end % i

H = hashtable;

for i=1:length(headers),
	H(headers{i}) = i;
end

% 2. M, C, R

M = zeros(lc-1, length(headers));
C = cell(lc-1, length(headers)); 
R = cell(lc-1,1);

i = 0;
while ~feof(fid),
	line = fgetl(fid);
	i = i + 1;

	tokens = tokenizer(line, delimiter);
	for j=1:length(tokens),
		M(i,j) = str2double(tokens{j});
		C{i,j} = tokens{j};
    end

	R{i} = line;
end

for j=1:length(tokens),
	C2{j} = {C{:,j}}'; 
end
C = C2;

fclose(fid);

% S

K = keys(H);
for i = 1:length(K),
	if isnan(M(1,i)),
		eval(sprintf('S.%s = C{i};', K{i}));
	else,
		eval(sprintf('S.%s = M(:,i);', K{i}));
	end
end % i
