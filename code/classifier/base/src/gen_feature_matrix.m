function [out] = gen_feature_matrix(varargin)
% Start Cache
[success, INTERN_cache_desc, varargin] = cache_enter(varargin);
if (success) == 1, out = INTERN_cache_desc; return; end;


data = load_cached_object(varargin{1});
cond_name = varargin{2};
bands = varargin{3};
features = varargin{4};
normalize = varargin{5};

% id_dict
cond = data.(cond_name);
id_dict.subjects = unique(data.subject);
id_dict.blocks = unique(data.block);
id_dict.conds = unique(cond);

data.id_dict = id_dict;

%%%% generate task level features %%%%

% H
H = {};
channels = {};
for i = 1:length(bands) - 1,
    channels{i} = sprintf('band_%d_%d', bands(i), bands(i+1)-1);
end
for c = 1:length(channels),
    channel = channels{c};
    for f = 1:length(features),
        feature = features{f};
        H{end + 1} = sprintf('%s_%s', feature, channel);
    end
end

H = [H 'LATENCY' 'AOA'];
FEAT_RANGE = 1:length(H);

H = [H 'SUBJECT' 'BLOCK' 'TASK' 'COND' 'INTERN_QUALITY' 'ACCEPT'];
TASK_RANGE = (max(FEAT_RANGE)+1):(length(H));

H = upper(H);
list = H; enum = 1; enum_list;

% M
M = NaN(length(data.start_time), length(H));
for t = 1:length(data.start_time),
    for c = 1:length(channels),
        feats = data.higher_order_features{t};
        if isempty(feats)
            continue;
        end
        fns = fieldnames(feats);
        for f = 1:length(fns)
            idx = strcmp(H, upper(fns{f}));
            M(t, idx) = feats.(fns{f});
        end
    end
    M(t, SUBJECT) = strlocate(id_dict.subjects, data.subject{t});
    M(t, BLOCK) = strlocate(id_dict.blocks, data.block{t});
    M(t, TASK) = data.task{t};
    if (iscell(cond)),
        M(t, COND) = strlocate(id_dict.conds, cond{t});
    else
        M(t, COND) = find(id_dict.conds == cond(t));
    end
    M(t, INTERN_QUALITY) = data.INTERN_QUALITY{t};
    M(t, LATENCY) = data.latency(t);
    M(t, AOA) = data.aoa(t);
end

% Normalize all data per subject
if normalize,
    subjects = unique(M(:, SUBJECT));
    for s = 1:length(subjects),
        I_subject = M(:, SUBJECT) == s;
        M_feat_data = M(I_subject, FEAT_RANGE);
        % compute the zscore ignoring nan
        M(I_subject, FEAT_RANGE) = (M_feat_data - repmat(nanmean(M_feat_data), size(M_feat_data, 1), 1)) ...
            ./ repmat(nanstd(M_feat_data), size(M_feat_data, 1), 1);
    end
end

data.feat_H = H(FEAT_RANGE);
data.feat_M = M(:, FEAT_RANGE);
data.task_H = H(TASK_RANGE);
data.task_M = M(:, TASK_RANGE);

%%%% generate epoch level features %%%%

% H
epoch_headers = {};
for c = 1:length(channels),
    channel = channels{c};
    epoch_headers{end + 1} = sprintf('%s', channel);
end

epoch_headers = upper(epoch_headers);

% M
epoch_M = NaN(0, length(epoch_headers));
epoch_idx = NaN(0, 1);
for t = 1:length(data.epochs),
    epochs = data.epochs{t};
    for e = 1:length(epochs),
        if ~isfield(epochs(e), 'features') || isempty(epochs(e).features)
            continue;
        end
        feats = epochs(e).features;
        epoch_M(end+1, :) = NaN(1, length(epoch_headers));
        epoch_idx(end+1, 1) = t;
        fns = fieldnames(feats);
        for f = 1:length(fns)
            idx = strcmp(epoch_headers, upper(fns{f}));
            epoch_M(end, idx) = feats.(fns{f});
        end
    end
end

data.epoch.H = epoch_headers;
data.epoch.M = epoch_M;
data.epoch.idx = epoch_idx;

out = data;

% End Cache
out = cache_exit(INTERN_cache_desc, out);
end
