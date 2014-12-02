
clear;

load('../../../data-lg/subjects.mat');

subjectids = subjects.subjectid;

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
