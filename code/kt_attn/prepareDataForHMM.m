
clear;

load('allSubjects.mat');

%% Build sequences

MIN_SEQUENCE_LENGTH = 5;
PERCENT_GOOD_EEG_THRESH = .25;

goodSubjects = {};
for sid=subjectids'
    fprintf('%s -- ', char(sid));
    filename = char(strcat('subjects/', sid, '.mat'));
    seqs_filename = char(strcat('subjects/', sid, '_sequences.mat'));
    
    load(filename, 'data');
    
    uniqueWords = unique(data.stim);
    numUniqueWords = length(uniqueWords);
    
    % build word sequences
    words = cell(numUniqueWords,1);
    total_encounters = 0;
    accept = cell(numUniqueWords,1);
    timeelapsed = cell(numUniqueWords,1);
    all_timeelapsed = [];
    attention = cell(numUniqueWords,1);
    all_attention = [];
    meditation = cell(numUniqueWords,1);
    hasfulleeg = zeros(numUniqueWords,1);
    seqsIdx = 1;
    for i=1:numUniqueWords
        word = uniqueWords{i};
        wordIdx = strcmp(data.stim, word);

        if sum(wordIdx)>MIN_SEQUENCE_LENGTH
            words{seqsIdx} = word;
            total_encounters = total_encounters + sum(wordIdx);
            asrResult = data.fluent(wordIdx);
            % ensure that asrResult vals are in (1,2)  
            if length(unique(asrResult))==3 || ...
                    (any(asrResult==0) && any(asrResult==2))
%                 fprintf('fluent contains both 0 and 2!!\n');
                continue;
            end
            if isempty(find(asrResult==2))
                asrResult = double(asrResult) + ones(size(asrResult));
            else
                a=1;
            end
            accept{seqsIdx} = asrResult;
            
            timeelapsed{seqsIdx} = data.timeelapsedms(wordIdx);
            all_timeelapsed = [all_timeelapsed; timeelapsed{seqsIdx}];
            attentionForSequence = data.attention(wordIdx);
            attention{seqsIdx} = attentionForSequence;
            all_attention = [all_attention; attentionForSequence(attentionForSequence~=-1)];
            meditation{seqsIdx} = data.meditation(wordIdx);
            if ~any(attention{seqsIdx}==-1)
                hasfulleeg(seqsIdx) = 1;
            end
            seqsIdx = seqsIdx + 1;
        end
    end
    
    words(seqsIdx:end) = [];
    accept(seqsIdx:end) = [];
    timeelapsed(seqsIdx:end) = [];
    attention(seqsIdx:end) = [];
    meditation(seqsIdx:end) = [];
    hasfulleeg(seqsIdx:end) = [];
    
    percentEncountersWithEEG = length(all_attention)/total_encounters;
    
    if seqsIdx == 1
        fprintf('NO sequences\n');
        continue;
    elseif percentEncountersWithEEG < PERCENT_GOOD_EEG_THRESH
        fprintf('Low EEG percentage\n');
        continue;        
    else
        fprintf('%d sequences (%d with full eeg)\n', seqsIdx, sum(hasfulleeg));
        goodSubjects{end+1} = sid;
    end

    sequences.subjectid = sid;
    sequences.words = words;
    sequences.accept = accept;
    sequences.timeelapsed = timeelapsed;
    sequences.meanTimeElapsed = mean(all_timeelapsed);
    sequences.varianceTimeElapsed = std(double(all_timeelapsed));
    sequences.attention = attention;
    sequences.meanAttention = mean(all_attention);
    sequences.varianceAttention = std(double(all_attention));
    sequences.percentEncountersWithEEG = percentEncountersWithEEG;
    sequences.meditation = meditation;
    sequences.hasfulleeg = hasfulleeg;
    
    save(seqs_filename, 'sequences');
    
end

%% Plot percentage of encounters with good EEG data for each subject

N = length(subjectids);
figure; hold on;
for i=1:N
    sid=subjectids{i};
    seqs_filename = char(strcat('subjects/', sid, '_sequences.mat'));
    load(seqs_filename);
    plot(i, sequences.percentEncountersWithEEG, 'o');
end
plot([1 N], [PERCENT_GOOD_EEG_THRESH PERCENT_GOOD_EEG_THRESH], 'k--');

%% Save good subjects

subjectids=goodSubjects;
save('subjects', 'subjectids');
