function [out] = apply_feature_selector(varargin)

data = varargin{1};
trained_selector = varargin{2};

feature_selector = trained_selector.method;
input = trained_selector.input;

list = data.task_H; enum = 1; enum_list;

switch feature_selector.algorithm,
    case 'pca',
        dimensions = feature_selector.dimensions;
        coef = pca(input.feat_M);
        coef = coef(:, 1:min(size(coef,2), dimensions));
        data.feat_M = data.feat_M * coef;
        
        feat_H = cell(1, dimensions);
        for i = 1:dimensions,
            feat_H{i} = sprintf('PCA_FEAT_%d', i);
        end
        data.feat_H = feat_H;
    case 'rank',
        selected = rankfeatures(transpose(input.feat_M), ...
            transpose(input.task_M(:, COND)), ...
            'NumberOfIndices', ...
            feature_selector.num_selected);
        data.feat_H = data.feat_H(selected);
        data.feat_M = data.feat_M(:, selected);
    otherwise,
        fprintf('select_features: error: feature selector %s is not supported\n', ...
            feature_selector.algorithm);
end

out = data;

end