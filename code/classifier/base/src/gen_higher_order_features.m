function [ out ] = gen_higher_order_features(varargin)
% Start Cache
[success, INTERN_cache_desc, varargin] = cache_enter(varargin);
if (success) == 1, out = INTERN_cache_desc; return; end;

data = varargin{1};
features = varargin{2};

for t = 1:length(data.start_time)
    if isempty(data.epochs{t}) || ~isfield(data.epochs{t}(1), 'features')
        data.higher_order_features{t} = [];
        continue;
    end
    epoch_features = data.epochs{t}(1).features;
    for e = 2:numel(data.epochs{t})
        epoch = data.epochs{t}(e);
        epoch_features(end+1) = epoch.features;
    end
    
    task_features = struct;
    fields = fieldnames(epoch_features);
    for i = 1:length(features),
        feat = features{i};
        for j = 1:length(fields),
            fn = fields{j};
            task_features.(sprintf('%s_%s', feat, fn)) = ...
                compute_feature(feat, [epoch_features.(fn)]');
        end
    end
    data.higher_order_features{t} = task_features;
end

out = data;

% End Cache
out = cache_exit(INTERN_cache_desc, out);
end

function value = compute_feature(feat, xs)
    switch feat,
        case 'mean',
            value = nanmean(xs, 1);
        case 'var',
            value = nanvar(xs,[],1);
        case 'min',
            value = nanmin(xs,[],1);
        case 'max',
            value = nanmax(xs,[],1);
        case 'sum',
            value = sum(~isnan(xs),1);
        case 'skew',
            value = skewness(xs,[],1);
        case 'kurtosis',
            value = kurtosis(xs,[],1);
    end
end