
clear;

load('subjects.mat');

MIN_SEQUENCE_LENGTH = 5;

for sid=subjectids'
    fprintf('%s -- ', char(sid));
    filename = char(strcat('subjects/', sid, '.mat'));
    seqs_filename = char(strcat('subjects/', sid, '_sequences.mat'));
    
    load(filename, 'data');
    
    uniqueWords = unique(data.stim);
    numUniqueWords = length(uniqueWords);
    
    words = cell(numUniqueWords,1);
    accept = cell(numUniqueWords,1);
    timeelapsed = cell(numUniqueWords,1);
    attention = cell(numUniqueWords,1);
    meditation = cell(numUniqueWords,1);
    hasfulleeg = zeros(numUniqueWords,1);
    seqsIdx = 1;
    for i=1:numUniqueWords
        word = uniqueWords{i};
        wordIdx = strcmp(data.stim, word);

        if sum(wordIdx)>MIN_SEQUENCE_LENGTH
            words{seqsIdx} = word;
            accept{seqsIdx} = data.fluent(wordIdx);
            timeelapsed{seqsIdx} = data.timeelapsedms(wordIdx);
            attention{seqsIdx} = data.attention(wordIdx);
            meditation{seqsIdx} = data.meditation(wordIdx);
            if ~any(attention{seqsIdx}==-1)
                hasfulleeg(seqsIdx) = 1;
            end
            seqsIdx = seqsIdx + 1;
        end
    end
    
    if seqsIdx == 1
        fprintf('NO sequences\n');
        continue;
    else
        fprintf('%d sequences (%d with full eeg)\n', seqsIdx, sum(hasfulleeg));
    end

    words(seqsIdx:end) = [];
    accept(seqsIdx:end) = [];
    timeelapsed(seqsIdx:end) = [];
    attention(seqsIdx:end) = [];
    meditation(seqsIdx:end) = [];
    hasfulleeg(seqsIdx:end) = [];
    
    sequences.subjectid = sid;
    sequences.words = words;
    sequences.accept = accept;
    sequences.timeelapsed = timeelapsed;
    sequences.attention = attention;
    sequences.meditation = meditation;
    sequences.hasfulleeg = hasfulleeg;
    
    save(seqs_filename, 'sequences');
    
end

load gong.mat;
soundsc(y);