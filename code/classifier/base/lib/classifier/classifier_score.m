function [acc, racc, eY, rank, confusion] = classifier_score(Y, S);

assertNaNScore = false;
assertScrambledClassIndices = false;

Y0 = Y;
S0 = S;

acc = NaN(length(Y),1);
racc = NaN(length(Y),1);
eY = NaN(length(Y),1);
rank = NaN(length(Y),1);
confusion = zeros(size(S,2), size(S,2));

% <assertion>
if sum(sum(isnan(S))) ~= 0,
	%warning('classifier_score.m:NaNScore', ...
		%'Scores have NaN values.');

	assertNaNScore = true;

	I_nan = sum(isnan(S),2) ~= 0;

	Y = Y0(~I_nan,:);
	S = S0(~I_nan,:);
end

if size(Y,1) < 1,
	%warning('classifier_score.m:emptyY', ...
		%'Empty Y');

	return;
end

if size(Y,1) ~= size(S,1) || length(unique(Y)) > size(S,2),
	warning('classifier_score.m:MismatchedYS', ...
		'Y and S are mismatched.');

	return;
end

if size(S,2) < 2,
	warning('classifier_score.m:NumberClassesLessThanTwo', ...
		'Number of distinct classes less than two.');

	return;
end

if length(unique(Y)) == size(S,2) && ~isequal(unique(Y)', 1:length(unique(Y))),
	warning('classifier_score.m:ScrambledClassIndices', ...
		'Class indices are scrambled.');

	return;

	%assertScrambledClassIndices = true;

	%uY = unique(Y);
	%for i = 1:length(Y),
		%Y2(i) = find(uY == Y(i));
	%end

	%Y = Y2';
end
% </assertion>

N = size(S,2);

[temp I] = sort(S * -1, 2);

eY = I(:,1);

acc = (Y == eY);

rank = zeros(length(Y), 1);
for i = 1:length(Y),
	rank(i) = find(I(i,:) == Y(i));
end

racc = 1 - (rank - 1) / (N - 1);

confusion = zeros(N);
for i = 1:length(Y),
	confusion(eY(i), Y(i)) = confusion(eY(i), Y(i)) + 1;
end

% <assertion>
if assertScrambledClassIndices,
	Y = Y0;

	uY = unique(Y);
	for i = 1:length(eY),
		eY2(i) = uY(eY(i));
	end

	eY = eY2';
end

if assertNaNScore,
	Y = Y0;
	S = S0;

	acc2 = NaN(length(Y),1);
	racc2 = NaN(length(Y),1);
	eY2 = NaN(length(Y),1);
	rank2 = NaN(length(Y),1);

	acc2(~I_nan,:) = acc;
	racc2(~I_nan,:) = racc;
	eY2(~I_nan,:) = eY;
	rank2(~I_nan,:) = rank;

	acc = acc2;
	racc = racc2;
	eY = eY2;
	rank = rank2;
end
% </assertion>
