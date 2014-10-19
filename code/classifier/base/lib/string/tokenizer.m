function [tokens] = tokenizer(str, delimiter)

k = 0;
tokens = {};

if nargin < 2,
   delimiter = '	 ';
end

while true
	[token, str] = strtok(str, delimiter);
	if isempty(token),  break;  end
	k = k + 1;
	tokens{k} = token;
end
