function [data,hmmData] = prepareData( task_file )


%% WARNING: this file is OLD!! Use prepareDataforHMM instead



%% WARNING: this file is OLD!! Use prepareDataforHMM instead



%% WARNING: this file is OLD!! Use prepareDataforHMM instead



%% WARNING: this file is OLD!! Use prepareDataforHMM instead



%% load input files

data=tdfread(task_file);

% convert chars to cell arrays
data.subject=cellstr(data.subject);
data.stim=cellstr(data.stim);

%% put it into form for HMM

MIN_SEQUENCE_LENGTH = 5;

hmmData = {};
subjectIDs = unique(data.subject);
for si = 1:length(subjectIDs)
    subjectID = subjectIDs{si};
    idx = strcmp(data.subject, subjectID);

    allWordsForSubject = data.stim(idx);
    % whether or not they got the word correct.
    % options here include: 'correct', 'latency', and 'fluent'
    allResultsForSubject_accept = data.correct(idx); 
    allResultsForSubject_latency = data.latency(idx); 
    allResultsForSubject_fluent = data.fluent(idx); 

    uniqueWords = unique(allWordsForSubject);
    
    % remove 'counts' at the end of some words
    uniqueWordsCleaned = {};
    for wi=1:length(uniqueWords)
        w = uniqueWords{wi};
        a = regexp(w,'\(\d+\)');
        if ~isempty(a)
            w = w(1:a-1);
        end
        uniqueWordsCleaned{wi} = w;
    end
    uniqueWords = unique(uniqueWordsCleaned);
    
    numUniqueWords = length(uniqueWords);

    words = cell(numUniqueWords,1);
    sequences_accept = cell(numUniqueWords,1);
    sequences_latency = cell(numUniqueWords,1);
    sequences_fluent = cell(numUniqueWords,1);

    seqsIdx = 1;
    for i=1:numUniqueWords
        word = uniqueWords{i};
        wordIdx = strcmp(allWordsForSubject, word);

        if sum(wordIdx)>MIN_SEQUENCE_LENGTH
            words{seqsIdx} = word;
            sequences_accept{seqsIdx} = allResultsForSubject_accept(wordIdx);
            sequences_latency{seqsIdx} = allResultsForSubject_latency(wordIdx);
            sequences_fluent{seqsIdx} = allResultsForSubject_fluent(wordIdx);
            seqsIdx = seqsIdx + 1;
        end
    end

    words(seqsIdx:end) = [];
    sequences_accept(seqsIdx:end) = [];
    sequences_latency(seqsIdx:end) = [];
    sequences_fluent(seqsIdx:end) = [];
    
    dataForSubject.subjectID = subjectID;
    dataForSubject.words = words;
    dataForSubject.accept = sequences_accept;
    dataForSubject.latency = sequences_latency;
    dataForSubject.fluent = sequences_fluent;
    hmmData{si} = dataForSubject;
end

end
