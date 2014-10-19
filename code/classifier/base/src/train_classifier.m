function out = train_classifier(varargin)

data = varargin{1};
classifier = varargin{2};

list = data.task_H; enum = 1; enum_list;
train.in = data.feat_M;
train.out = data.task_M(:, COND);
ok = 1;

if strcmp(classifier, 'svm'),
    options = optimset('maxiter',10000);
    try,
        c = svmtrain(train.in, train.out,...
            'kernel_function','linear','options', options, 'boxconstraint', 2);
    catch err,
        err
        ok = 0;
        if (strcmp(err.identifier, 'stats:svmtrain:NoConvergence')),
        else,
            error('run_classification.m\between\svm');
        end
    end
else
    c = trainClassifier(train.in, train.out, classifier);
end

c_struct.type = classifier;
c_struct.ok = ok;
if c_struct.ok == 1,
    c_struct.classifier = c;
end

out = c_struct;

end