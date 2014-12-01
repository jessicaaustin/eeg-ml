
clear;

load('subjects.mat', 'subjectids');

% sid='mjs6-8-2001-09-24'
for sid=subjectids'
    tic
    
    fprintf('%s -- ', char(sid));
    idx_filename = char(strcat('subjects/', sid, '_idx.mat'));
    rawdata_filename = char(strcat('../../../data-lg/userdata/', sid, '_aligned.xls'));
    filename = char(strcat('subjects/', sid, '.mat'));
    
    try
        rawdata = tdfread(rawdata_filename);
    catch err
        fprintf('ERR\n');
        continue;
    end
    load(idx_filename, 'idx');
    goodeegidx = zeros(size(rawdata.haseeg));
    goodeegidx(idx)=1;
    
    data.subjectid = sid;
    data.userid = rawdata.User_ID(idx(1));
    data.stim = cellstr(rawdata.stim);
    data.timeelapsedms = int32(rawdata.timeelapsedms);
    data.fluent = int8(rawdata.fluent);
    data.haseeg = goodeegidx;
    data.attention = int32(rawdata.Attention);
    data.attention(~goodeegidx) = -1;
    data.meditation = int32(rawdata.Meditation);
    data.meditation(~goodeegidx) = -1;
    
    save(filename, 'data');
    fprintf('%s\n', filename);
    
    toc
end

load gong.mat;
soundsc(y);