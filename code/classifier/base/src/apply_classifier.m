function out = apply_classifier(varargin)

data = varargin{1};
c = varargin{2};

feats = data.feat_M;
results = zeros(size(data.task_M, 1), length(data.id_dict.conds));

if strcmp(c.type, 'svm'),
    predictions = svmclassify(c.classifier, feats);
    for p = 1:numel(predictions),
        results(p, predictions(p)) = 1;
    end
else
    results = applyClassifier(feats, c.classifier);
end

out = results;

end

