function out = gen_cv_splits(varargin)
% Start Cache
[success, INTERN_cache_desc, varargin] = cache_enter(varargin);
if (success) == 1, out = INTERN_cache_desc; return; end;

data = varargin{1};
cv = varargin{2};
cv_subjects = varargin{3};

global DEBUG VERBOSE;

if exist('VERBOSE', 'var') && ~isempty(strfind(VERBOSE, '-printTicToc')),
	fprintf('FUNCTION %s, TOC %.2f min\n', mfilename, toc / 60);
end

H = data.task_H;
M = data.task_M;
list = H; enum = 1; enum_list;

splits = struct('train', {}, 'test', {});
subjects = unique(M(:, SUBJECT));
classes = length(data.id_dict.conds);

subject_splits = {};
switch cv_subjects,
    case 'within',
        for s = 1:length(subjects),
            subject_splits{end+1} = M(:,SUBJECT) == s;
        end
    case 'between',
        for s = 1:length(subjects),
            I = M(:,SUBJECT) ~= s;
            splits(end+1).train = find(I);
            splits(end).test = find(~I);
        end
    otherwise,
        subject_splits{end+1} = logical(ones(size(M, 1), 1));
end

if ~strcmp(cv_subjects, 'between')
    for s = 1:length(subject_splits)
        I = subject_splits{s};
        is = find(I);
        tasks = unique(M(I, TASK));
        switch cv,
            case 'leave-one-out',
                for i = 1:length(tasks),
                    t = tasks(i);
                    I_t = M(I, TASK) == t;
                    train = I;
                    train(is(I_t)) = 0;
                    test = zeros(length(I), 1);
                    test(is(I_t)) = 1;
                    splits(end+1).train = find(train);
                    splits(end).test = find(test);
                end
        end
    end
end

% remove underspecified cv folds
i = 1;
while i <= length(splits)
    cv = splits(i);
    if length(unique(M(cv.train, COND))) < classes,
        splits(i) = [];
        continue;
    end
    i = i + 1;
end

out = splits;

% End Cache
out = cache_exit(INTERN_cache_desc, out);
end

