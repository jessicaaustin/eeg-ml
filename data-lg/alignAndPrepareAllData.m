
clear;

% % load subjects
subjects=tdfread('subjectToUserId.xls');
subjects.subjectid=cellstr(subjects.subjectid);
save('subjects', 'subjects');

for sid=1:length(subjects.subjectid)
    
    subjectid=subjects.subjectid{sid};
    userid=subjects.userid(sid);
    
    disp(subjectid);
    
    asrfilename = char(strcat('userdata/', subjectid, '_asr.xls'));
    eegfilename = char(strcat('userdata/', subjectid, '_eeg.xls'));
    resultsfilename = char(strcat('userdata/', subjectid, '_aligned.xls'));
    
    % load ASR data
    asrfile=tdfread(asrfilename);
    N = size(asrfile.subject,1);
    for f=fieldnames(asrfile)'
        field = asrfile.(char(f));
        if ischar(field)
            asrfile.(char(f)) = cellstr(field);
        end
    end
    
    % load EEG data
    eegfile=tdfread(eegfilename);
    eegfields = fieldnames(eegfile)';
    for f=eegfields
        field = eegfile.(char(f));
        if ischar(field)
            eegfile.(char(f)) = cellstr(field);
        end
    end
    
    % create empty result set
    results = asrfile;
    for f=eegfields
        if iscell(eegfile.(char(f)))
            results.(char(f)) = repmat({' '}, N, 1);
        else
            results.(char(f)) = zeros(N,1);
        end
        results.timeelapsedms = zeros(N,1);
        results.haseeg = zeros(N,1);
    end
    
    % align data
    starttimes=datenum(asrfile.start_time);
    endtimes=datenum(asrfile.end_time);
    eegtimes=datenum(eegfile.Start_Time);
    for i=1:N
        t1 = starttimes(i);
        t2 = endtimes(i);
        idx = find(eegtimes>t1 & eegtimes<t2);
        if ~isempty(idx)
            for f=eegfields
                if iscell(eegfile.(char(f)))
                    results.(char(f)){i} = eegfile.(char(f)){idx(1)};
                else
                    results.(char(f))(i) = eegfile.(char(f))(idx(1));
                end
            end
            results.haseeg(i) = 1;
        end
        
        % store time elapsed
        results.timeelapsedms(i) = round((t2-t1)*24*60*60*1000);
        
        % clean up STIM
        stim = asrfile.stim{i};
        a = regexp(stim,'\(\d+\)');
        if ~isempty(a)
            stim = stim(1:a-1);
        end
        results.stim{i} = stim;
    end
    
    % write results
    for f=fieldnames(results)'
        field = results.(char(f));
        if iscell(field)
            results.(char(f)) = char(field);
        end
    end
    tdfwrite(resultsfilename, results);
end

load gong.mat;
soundsc(y);