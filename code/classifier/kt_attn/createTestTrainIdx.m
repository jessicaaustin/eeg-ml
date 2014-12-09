clear;
load('subjects.mat');
N = length(subjectids);

testTrainIdx = {};

for i=1:N
    subjectid=subjectids{i};
    fprintf('%s:\n', char(subjectid));

    % load sequences
    seqs_filename = char(strcat('subjects/', subjectid, '_sequences.mat'));
    load(seqs_filename);
    
    % split into testing and training set
    sN = length(sequences.accept);
    Ntrain = round(.9*sN);
    allIdx = (1:sN)';
    idx = randperm(sN);
    idxTrain = sort(allIdx(idx(1:Ntrain)));
    idxTest = sort(allIdx(idx(Ntrain+1:end)));
    testTrainIdx{i}.idxTrain = idxTrain;
    testTrainIdx{i}.idxTest = idxTest;
    
end

save('testTrainIdx', 'testTrainIdx');

