function [out] = read_task(varargin)

% Input:
%
% expt.task_file
% expt.cache
% expt.cache_task
%
% Output:
%
% H
% M
% expt.machines
% expt.subjects
% expt.blocks

global DEBUG VERBOSE;

if exist('VERBOSE', 'var') && ~isempty(strfind(VERBOSE, '-printTicToc')),
    fprintf('FUNCTION %s, TOC %.2f min\n', mfilename, toc / 60);
end

% Start Cache
[success, INTERN_cache_desc, varargin] = cache_enter(varargin);
if (success) == 1, out = INTERN_cache_desc; return; end;

task_file = varargin{1};

% load task
[M H C R task] = tablescan(task_file);
H = upper(keys(H));
list = H; enum = 1; enum_list;

fields = fieldnames(task(1));
fields = fields(...
    cellfun(@(f) ~(strcmp(f, 'start_time') | strcmp(f, 'end_time')), ...
    fields));
% TODO separate by cell/mat
cell_idx = cellfun(@(f) iscell(task.(f)), fields);
cell_fields = fields(cell_idx);
mat_fields = fields(~cell_idx);

idx = 1;
for t = 1:size(M, 1),
    task2.start_time{idx} = datenum(task.start_time{t});
    task2.end_time{idx} = datenum(task.end_time{t});
    task2.task{idx} = t;
    for i = 1:length(cell_fields)
        f = cell_fields{i};
        task2.(f){idx} = task.(f){t};
    end
    for i = 1:length(mat_fields)
        f = mat_fields{i};
        task2.(f)(idx) = task.(f)(t);
    end
    idx = idx + 1;
end

if ~iscell(task.subject)
    task2.subject = cellstr(num2str(task.subject)); % support non-numerical user_ids
end

out = task2;

% End Cache
out = cache_exit(INTERN_cache_desc, out);
end