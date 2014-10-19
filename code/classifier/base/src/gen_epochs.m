function [out] = gen_epochs(varargin)
% Start Cache
[success, INTERN_cache_desc, varargin] = cache_enter(varargin);
if (success) == 1, out = INTERN_cache_desc; return; end;

data = varargin{1};
epochs = varargin{2};

overlap = datenum(0,0,0,0,0,epochs.overlap);
seg_size = datenum(0,0,0,0,0,epochs.length);
min_size = datenum(0,0,0,0,0,epochs.min);

for t = 1:length(data.start_time)
    % GENERATE OSCILLATION BANDS
    task_epochs = struct('time', [], 'signal', []);
    idx = 1;
    cur_time = data.start_time{t};
    while cur_time < data.end_time{t} - min_size
        start_time = cur_time;
        end_time = start_time + seg_size;
        
        I = find(and(start_time <= data.RAWWAVE(t).time, data.RAWWAVE(t).time <= end_time));
        ts = data.RAWWAVE(t).time(I);
        xs = data.RAWWAVE(t).signal(I);
        
        if (length(ts) > 2)
            task_epochs(idx).time = ts;
            task_epochs(idx).signal = xs;
            idx = idx + 1;
        end
        
        cur_time = start_time + seg_size - overlap;
    end
    data.epochs{t} = task_epochs;
end

out = data;

% End Cache
out = cache_exit(INTERN_cache_desc, out);
end
