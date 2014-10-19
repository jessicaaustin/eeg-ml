function [data, results] = run_experiment( expt )
%RUN_EXPERIMENT Summary of this function goes here
%   Detailed explanation goes here

run_setup;

data = run_prepare_data(expt);

% cross-validated classification
cv_splits = gen_cv_splits(data, expt.cv, expt.cv_subjects);
cv_results = run_all_classification(data, cv_splits, expt.feature_selector, expt.classifier, expt.balance, 0);

% evaluate results
data = load_cached_object(data);
list = data.task_H; enum = 1; enum_list;
gold = data.task_M(:, COND);

results = aggregate_results(data, cv_splits, cv_results);
[n accuracy p] = evaluate_results(gold, results);

end

function results = aggregate_results(data, cv_splits, cv_results)
    cv_splits = load_cached_object(cv_splits);
    cv_results = load_cached_object(cv_results);
    results = zeros(size(data.task_M, 1), length(data.id_dict.conds));
    for i = 1:length(cv_splits)
        cv = cv_splits(i);
        results(cv.test, :) = cv_results{i};
    end
end


