
clear;

load('../../../data-lg/subjects.mat');

subjectids = subjects.subjectid;

%%

for sid=subjectids'
    fprintf('%s -- ', char(sid));
    filename = char(strcat('subjects/', sid, '_idx.mat'));
    
    try
        rawdata = tdfread(char(strcat('../../../data-lg/userdata/', sid, '_aligned.xls')));
    catch err
        continue;
    end
    signal = rawdata.PoorSignal;
    if ~isnumeric(signal)
        signal = str2double(cellstr(signal));
    end
    idx = find(signal<=100);
    
    if ~isempty(idx)
        fprintf('%d vals\n', length(idx));
        save(filename, 'idx');
    else
        fprintf('No data :(\n');
    end
    
end


%%

total_num_encounters = 0;
total_encounters_with_eeg = 0;
total_subjects_with_eeg = 0;

for sid=subjectids'
    filename = char(strcat('subjects/', sid, '_idx.mat'));
    if exist(filename, 'file')
        load(filename, 'idx');
        total_subjects_with_eeg = total_subjects_with_eeg + 1;
        total_num_encounters = total_num_encounters + max(idx);
        total_encounters_with_eeg = total_encounters_with_eeg + length(idx);
    end
end

total_subjects_with_eeg
total_encounters_with_eeg/total_num_encounters

    