function [scoreAcrossTrials] = repl2(Y, X, measure)

if nargin < 3,
	measure = 'correlation';
end

if isempty(Y),
	scoreAcrossTrials = NaN(1, size(X,2));
	return;
end

F = zeros(length(Y), 1);
uY = unique(Y);
for y = 1:length(uY),
	I = find(Y == uY(y));
	F(I) = 1:length(I);
end

nuY = length(unique(Y));
nuF = length(unique(F));

% sort by F x Y
%[temp I] = sort(F * nuY * 2 + Y); % bug
[temp I] = sort(F * max(Y) + Y);
	X = X(I,:);
	Y = Y(I);
	F = F(I);

for i = 1:size(X,2);
	M = reshape(X(:,i), nuY, nuF)';
	scoreAcrossTrials(i) = mean(1 - pdist(M, measure));
	%reshape(Y, nuY, nuF)'
	%keyboard;
end
