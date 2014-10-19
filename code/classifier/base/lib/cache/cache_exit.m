function out = cache_exit(INTERN_cache_desc, out)
    global CACHE;
    if exist('CACHE', 'var') && ~strcmp(CACHE, ''),
        write_cache(INTERN_cache_desc, out);
        out = INTERN_cache_desc;
    end
end

function write_cache(desc, data)
    global CACHE INDEX;
    
    index_loc = sprintf('%s/index.mat', CACHE);
    if isempty(INDEX),
        if exist(index_loc, 'file'),
            load(index_loc);
            INDEX = index;
        else
            INDEX = {};
        end
    end
    
    INDEX{end+1} = desc;
    save(sprintf('%s/%d.mat', CACHE, length(INDEX)), 'data', '-v7.3');
    index = INDEX;
    save(sprintf('%s/index.mat', CACHE), 'index', '-v7.3');
end
