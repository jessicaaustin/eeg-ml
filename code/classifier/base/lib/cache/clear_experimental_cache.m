function clear_experimental_cache()
    global CACHE INDEX;
    
    id = -1;
    
    if isempty(INDEX),
        if exist(sprintf('%s/index.mat', CACHE), 'file'),
            load(sprintf('%s/index.mat', CACHE));
            INDEX = index;
        else
            return;
        end
    end
    
    i = 1;
    while i <= numel(index)
        desc = index{i};
        if is_experimental(desc) == 0
            i = i + 1;
            continue;
        end
        if numel(index) ~= i
            movefile(sprintf('%s/%d.mat', CACHE, numel(index)),...
                sprintf('%s/%d.mat', CACHE, i));
            index(i) = index(end);
        end
        index(end) = [];
    end
    INDEX = index;
end

function experimental = is_experimental(desc)
    experimental = 0;
    if isfield(desc, 'INTERN_experimental')
        experimental = 1;
        return;
    end
    if isfield(desc, 'vars')
        for i = 1:numel(desc.vars)
            experimental = experimental + is_experimental(desc.vars{i});
        end
    end
end
        
            
        