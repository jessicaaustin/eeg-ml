function [out] = run_all_classification(varargin)
% Start Cache
[success, INTERN_cache_desc, varargin] = cache_enter(varargin);
if (success) == 1, out = INTERN_cache_desc; return; end;

data = load_cached_object(varargin{1});
cv_splits = load_cached_object(varargin{2});
feature_selector = varargin{3};
classifier = varargin{4};
balance = varargin{5};
epoch_based = varargin{6};

results = cell(length(cv_splits), 1);
for i = 1:length(cv_splits),
    cv = cv_splits(i);
    
    cv.train = balance_data(cv.train, data, balance);
    if (epoch_based),
        cv_data = gen_cv_epoch_data(cv, data);
    else,
        cv_data = gen_cv_data(cv, data);
    end

    [cv_data.train, trained_selector] = train_feature_selector(cv_data.train, feature_selector);
    cv_data.test = apply_feature_selector(cv_data.test, trained_selector);
    
    trained_classifier = train_classifier(cv_data.train, classifier);
    if trained_classifier.ok,
        cv_results = apply_classifier(cv_data.test, trained_classifier);
    else,
        cv_results = ones(size(cv_data.test.task_M, 1), length(data.id_dict.conds)) * 0.5;
    end
    
    if (epoch_based),
        cv_results = epoch_vote(cv, data, cv_results);
    end
    results{i} = cv_results;
end

out = results;

% End Cache
out = cache_exit(INTERN_cache_desc, out);
end

function [out] = epoch_vote(cv, data, cv_results),
ecv.test = find(ismember(data.epoch.idx, cv.test));
t.ecv.test = data.epoch.idx(ecv.test);
voted_results = zeros(length(cv.test), size(cv_results, 2));
for i = 1:length(cv.test),
    I = (t.ecv.test == cv.test(i));
    results_i = cv_results(I, :);
    % TODO: make this work for more than 2 classes
    class2 = sum(results_i(:, 1) < results_i(:, 2));
    voted_results(i, 1 + (class2 > (size(results_i, 1) / 2))) = 1;
end
out = voted_results;
end

function [out] = gen_cv_epoch_data(varargin)
cv = varargin{1};
data = varargin{2};

% fill features
ecv.train = find(ismember(data.epoch.idx, cv.train));
ecv.test = find(ismember(data.epoch.idx, cv.test));
cv_data.train.feat_M = data.epoch.M(ecv.train, :);
cv_data.test.feat_M = data.epoch.M(ecv.test, :);
cv_data.train.feat_H = data.epoch.H;
cv_data.test.feat_H = data.epoch.H;

% fill task info
t.ecv.train = data.epoch.idx(ecv.train);
t.ecv.test = data.epoch.idx(ecv.test);
cv_data.train.task_M = data.task_M(t.ecv.train, :);
cv_data.test.task_M = data.task_M(t.ecv.test, :);
cv_data.train.task_H = data.task_H;
cv_data.test.task_H = data.task_H;
cv_data.train.id_dict = data.id_dict;
cv_data.test.id_dict = data.id_dict;
out = cv_data;

end

function [out] = gen_cv_data(varargin)
cv = varargin{1};
data = varargin{2};

% get data
cv_data.train.feat_M = data.feat_M(cv.train, :);
cv_data.test.feat_M = data.feat_M(cv.test, :);
cv_data.train.feat_H = data.feat_H;
cv_data.test.feat_H = data.feat_H;
cv_data.train.task_M = data.task_M(cv.train, :);
cv_data.test.task_M = data.task_M(cv.test, :);
cv_data.train.task_H = data.task_H;
cv_data.test.task_H = data.task_H;
cv_data.train.id_dict = data.id_dict;
cv_data.test.id_dict = data.id_dict;
out = cv_data;

end