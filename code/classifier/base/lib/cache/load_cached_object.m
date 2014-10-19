function data = load_cached_object(desc)
    if (~isfield(desc, 'INTERN_name')),
        data = desc;
        return; % this is not a cache descriptor
    end
    
    global CACHE;
    id = find_object(desc);
    if id == -1,
        err('object should be cached but is not found');
    end
    fname = sprintf('%s/%d.mat', CACHE, id);
    if exist(fname, 'file'),
        load(fname);
    end
end

function id = find_object(desc)
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
    
    id = find(cellfun(@(s) isequal(s, desc), INDEX));
    if ~isempty(id)
        return;
    end
end