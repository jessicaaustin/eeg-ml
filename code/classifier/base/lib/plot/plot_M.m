function h = plot_M(M, varargin)

% ------------------------------------------------
% Options

for i = 1:size(M,2), XLabels{i} = 'X'; end;
for i = 1:size(M,1), YLabels{i} = 'Y'; end;
Title = '';

angle = 90;
showSimilarity = 0;
showDistinctivenessCorrelation = 0;

saveFig = 0;
printFig = 0;

sortByRow = 0;
sortByColumn = 0;

saveFig = 0;
printFig = 0;

args = varargin;
for i = 1:2:length(args),
	eval(sprintf('%s = args{i+1};', args{i}));
end

% ------------------------------------------------
% Special processing

if showSimilarity,
	M = 1 - squareform(pdist(M, 'euclidean'));
	XLabels = YLabels;
end

if showDistinctivenessCorrelation,
	Xdistinctiveness = 1 ./ sum(M);
	Xcorrelation = mean(1 - squareform(pdist(M', 'correlation')));
	for i = 1:length(XLabels),
		XLabels{i} = sprintf('%s (%.3f, %.3f)', XLabels{i}, Xdistinctiveness(i), Xcorrelation(i));
	end
end

XLabels = escape(XLabels);
YLabels = escape(YLabels);

% ------------------------------------------------
% Sort

if sortByRow,
	[Y I] = sort(M, 2, 'descend');
	M = M(:,I(sortByRow,:));
	XLabels = XLabels(I(sortByRow,:));

elseif sortByColumn,
	[Y I] = sort(M, 1, 'descend');
	M = M(I(:,sortByColumn),:);
	YLabels = YLabels(I(:,sortByColumn));
end

% ------------------------------------------------
% plot the heat diagram

%M

%h = figure;
%h = clf;
h = gcf;
datacursormode on;

if exist('clims', 'var'),
	imagesc(M, clims);
else,
	imagesc(M);
end

%colorbar;
set(gca, 'XTick',      1:length(XLabels));
set(gca, 'XTickLabel', numberStrings(XLabels));
set(gca, 'YTick',      1:length(YLabels));
set(gca, 'YTickLabel', numberStrings(YLabels));

title(escape(Title));

%if exist('angle', 'var') && size(M,2) > 1,
if size(M,2) > 1,
	%xticklabel_rotate([], angle);
	rotateticklabel(gca,angle);
end

if saveFig,
	saveas(h, sprintf('%s.fig', Title));
end

if printFig,
	print_ccbi(h);
end

return;

% ------------------------------------------------
% Events

dcm_obj = datacursormode(h);
set(dcm_obj, 'UpdateFcn', {@figure1_UpdateFcn, XLabels, YLabels});

hSortByRow_Edit = uicontrol(gcf, ...
		'Style', 'Edit', ...
		'String', num2str(sortByRow), ...
		'Position', [0,0,30,30], ...
		'Callback',{@sortByRow_Callback, M, XLabels, YLabels, Title});

hSortByColumn_Edit = uicontrol(gcf, ...
		'Style', 'Edit', ...
		'String', num2str(sortByColumn), ...
		'Position', [31,0,30,30], ...
		'Callback',{@sortByColumn_Callback, M, XLabels, YLabels, Title});

hTranspose_Button = uicontrol(gcf, ...
		'Style', 'PushButton', ...
		'String', 'Transpose', ...
		'Position', [61,0,90,30], ...
		'Callback',{@transpose_Callback, M, XLabels, YLabels, Title});

function txt = figure1_UpdateFcn(source, event, XLabels, YLabels)
	pos = get(event,'Position');
	x = pos(2);
	y = pos(1);
	M = get(get(event, 'Target'), 'CData');

	txt = sprintf('M(%d, %d) = %2.4f\n%s\n%s', x, y, M(x,y), YLabels{x}, XLabels{y});

function sortByRow_Callback(source, event, M, XLabels, YLabels, Title)
	plot_M(M, 'XLabels', XLabels, 'YLabels', YLabels, 'Title', Title, ...
			'sortByRow', str2num(get(source, 'String')));

function sortByColumn_Callback(source, event, M, XLabels, YLabels, Title)
	plot_M(M, 'XLabels', XLabels, 'YLabels', YLabels, 'Title', Title, ...
			'sortByColumn', str2num(get(source, 'String')));

function transpose_Callback(source, event, M, XLabels, YLabels, Title)
	plot_M(M', 'XLabels', YLabels, 'YLabels', XLabels, 'Title', Title);
