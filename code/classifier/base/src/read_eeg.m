function [out] = read_eeg(varargin)

% Input:
%
% expt.eeg_files
%

global DEBUG VERBOSE;

if exist('VERBOSE', 'var') && ~isempty(strfind(VERBOSE, '-printTicToc')),
	fprintf('FUNCTION %s, TOC %.2f min\n', mfilename, toc / 60);
end

% Start Cache
[success, INTERN_cache_desc, varargin] = cache_enter(varargin);
if (success) == 1, out = INTERN_cache_desc; return; end;

eeg_file = varargin{1};
sigqual_sample = varargin{2};

%% Load
[M, H, C, R, task] = tablescan(eeg_file);

H_SIGNAL = upper(keys(H));
list = H_SIGNAL; enum = 1; enum_list;

% H_SIGNAL, M_SIGNAL, M_RAW

% Create the matrix M_SIGNAL, which contains all data from the EEG file except raw signal in numeric form
% (Machine name replaced with index of machine in list of machines, dates replaced with date strings, etc.)
% Form of M_SIGNAL: first column has subject ID, second/third have start/end time, remaining ones have EEG feature data
% This version of the loop requires the Excel file to contain the numeric date with just : and / as separators
% (no extra letters or AM/PM), but runs faster than the other one

M0_SIGNAL = NaN(size(M,1), 4);
M0_RAW = cell(size(M,1), 1);
if ~iscell(task.subject)
    task.subject = cellstr(num2str(task.subject)); % support non-numerical user_ids
end
subjects = unique(task.subject);
for t = 1:size(M,1),
	M0_RAW{t} = str2num(task.rawwave{t});

    M0_SIGNAL(t, START_TIME) = datenum(task.start_time{t});
    M0_SIGNAL(t, END_TIME) = datenum(task.end_time{t});
    M0_SIGNAL(t, SUBJECT) = strlocate(subjects, task.subject{t});
	M0_SIGNAL(t, SIGQUAL) = M(t, SIGQUAL);
end % t

all_data = {};
for s = 1:length(subjects),
    % grab this subject's data
    I = M0_SIGNAL(:,SUBJECT) == s;
	M_SIGNAL = M0_SIGNAL(I,:);
	M_RAW = M0_RAW(I);
    
    % Convert rawwave matrix format into linear stream
    M_RAW_SIGNAL = cell2mat(M_RAW')';
    M_RAW_TIME = nan(size(M_RAW_SIGNAL));
    idx = 1;
    for t = 1:size(M_RAW, 1),
        nsample = length(M_RAW{t});
        if M0_SIGNAL(t, SIGQUAL) > sigqual_sample,
            % nan out bad samples
            M_RAW_SIGNAL(idx:idx + nsample - 1) = NaN;
        else
            % if the sample is good, add time data
            T = M_SIGNAL(t,START_TIME) + (0:nsample)' .* ( M_SIGNAL(t,END_TIME) - M_SIGNAL(t,START_TIME) ) / nsample;
            M_RAW_TIME(idx:idx + nsample - 1) = T(1:end-1);
        end
        idx = idx + nsample;
    end % t
    
    % remove all nans (bad values)
    M_RAW_SIGNAL = M_RAW_SIGNAL(isfinite(M_RAW_SIGNAL(:, 1)), :);
    M_RAW_TIME = M_RAW_TIME(isfinite(M_RAW_TIME(:, 1)), :);
    
    % write data to the struct
    data.subject = subjects{s};
    data.rawwave.time = M_RAW_TIME;
    data.rawwave.signal = M_RAW_SIGNAL;
    
    all_data{end + 1} = data;
end

out = all_data;

% End Cache
out = cache_exit(INTERN_cache_desc, out);
end
