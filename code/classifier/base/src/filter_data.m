function out = filter_data(varargin)
% Start Cache
[success, INTERN_cache_desc, varargin] = cache_enter(varargin);
if (success) == 1, out = INTERN_cache_desc; return; end;


data = varargin{1};
sigqual_percentage = varargin{2};

global DEBUG VERBOSE;

if exist('VERBOSE', 'var') && ~isempty(strfind(VERBOSE, '-printTicToc')),
	fprintf('FUNCTION %s, TOC %.2f min\n', mfilename, toc / 60);
end

list = data.task_H; enum = 1; enum_list;

M = data.task_M;
% filter by signal
I = M(:,INTERN_QUALITY) < sigqual_percentage;
M(I,:) = NaN;
filter.sigqual = [size(M(~isnan(M(:, COND)),:)) nanfreq(M(:, COND))];

% filter out rows with NaN features
I = any(isnan(data.feat_M), 2) | I;
M(I,:) = NaN;
filter.nan = [size(M(~isnan(M(:, COND)),:)) nanfreq(M(:, COND))];

filter

data.task_M = data.task_M(~I, :);
data.feat_M = data.feat_M(~I, :);

out = data;

% End Cache
out = cache_exit(INTERN_cache_desc, out);
end
