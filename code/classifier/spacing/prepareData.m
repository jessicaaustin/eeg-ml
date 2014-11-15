function [data,words,sequences] = prepareData( expt )
%% Experiment

% load input files
data = read_task(expt.task_file);
data = load_cached_object(data);

% put it into form for HMM

% TODO: load for all subjects
subjectID = 's2';
idx = strcmp(data.subject, subjectID);

allWordsForSubject = data.stim(idx);
allASRresultForSubject = data.accept(idx);

uniqueWords = unique(allWordsForSubject);
numUniqueWords = length(uniqueWords);

words = cell(numUniqueWords);
sequences = cell(numUniqueWords);

seqsIdx = 1;
for i=1:numUniqueWords
    word = uniqueWords{i};
    wordIdx = strcmp(allWordsForSubject, word);
    
    if sum(wordIdx)>5
        seq = allASRresultForSubject(wordIdx);

        words{seqsIdx} = word;
        sequences{seqsIdx} = seq;
        seqsIdx = seqsIdx + 1;
    end
end

words(seqsIdx:end) = [];
sequences(seqsIdx:end) = [];

end
