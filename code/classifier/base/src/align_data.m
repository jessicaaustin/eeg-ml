function [out] = align_data(varargin)

% Start Cache
[success, INTERN_cache_desc, varargin] = cache_enter(varargin);
if (success) == 1, out = INTERN_cache_desc; return; end;

task_data = varargin{1};
eeg_data = varargin{2};

global DEBUG VERBOSE;

if exist('VERBOSE', 'var') && ~isempty(strfind(VERBOSE, '-printTicToc')),
    fprintf('FUNCTION %s, TOC %.2f min\n', mfilename, toc / 60);
end

%% Create M
eeg = eeg_data{1}; % load something

% init data arrays
task_data.('INTERN_QUALITY') = num2cell(zeros(length(task_data.start_time), 1));
aligned_rawwave.time = cell(length(task_data.start_time), 1);
aligned_rawwave.signal = cell(length(task_data.start_time), 1);

for t = 1:length(task_data.start_time),
    subject = task_data.subject{t};
    start_time = task_data.start_time{t};
    end_time = task_data.end_time{t};
    
    % get EEG data
    if ~strcmp(eeg.subject, subject),
        found = 0;
        for i = 1:length(eeg_data),
            e = eeg_data{i};
            if strcmp(e.subject, subject),
                eeg = e;
                found = 1;
            end
        end
        if ~found,
            continue;
        end
    end
    
    % align rawwave
    I = find(and(start_time <= eeg.rawwave.time, eeg.rawwave.time <= end_time));
    ts = eeg.rawwave.time(I);
    xs = eeg.rawwave.signal(I);
    task_data.('INTERN_QUALITY'){t} = length(xs); % TODO change this
    
    % store rawwave
    aligned_rawwave.time{t} = ts;
    aligned_rawwave.signal{t} = xs;
end

task_data.('RAWWAVE') = struct('time', aligned_rawwave.time, 'signal', aligned_rawwave.signal); 

out = task_data;

% End Cache
out = cache_exit(INTERN_cache_desc, out);
end

