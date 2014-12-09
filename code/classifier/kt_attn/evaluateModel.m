clear;
% close all;

addpath('../../common');
addpath('../../HMM_mat');
addpath('../../HMM_mat_ext');

load('subjects.mat');
N = length(subjectids);
load('testTrainIdx');

rng('default');
rng(3);

%% Estimate Params

allActualAsrObservation = [];
allEstAsrObservation_KT = [];
allEstAsrObservation_KTAttn = [];

for i=1:N
    subjectid=subjectids{i};
    fprintf('%s:\n', char(subjectid));

    % load sequences
    seqs_filename = char(strcat('subjects/', subjectid, '_sequences.mat'));
    load(seqs_filename);
    
    % split into testing and training set
    idxTrain = testTrainIdx{i}.idxTrain;
    idxTest = testTrainIdx{i}.idxTest;

    % estimate params on training set
    [~,p_KT,A_KT,B_KT] = estimateParamsForSubject(sequences, idxTrain, 'KT');
    [~,p_KTAttn,A_KTAttn,B_KTAttn,C,D] = estimateParamsForSubject(sequences, idxTrain, 'KTAttn');
    
    % generate observations on testing set
    for si = idxTest(:)'
        
        actualAsrObservations = sequences.accept{si};
        actualEEGObservations = thresholdAndFillAttention(sequences.attention{si}, sequences);
        v0 = actualAsrObservations(1);
        w0 = actualEEGObservations(1);
        T = length(actualAsrObservations);
        
        allActualAsrObservation = [allActualAsrObservation;
                                   actualAsrObservations];
        
        % estimate observations
        asrObservations = runModelForward(p_KT,A_KT,B_KT,[],[],v0,[],T, 'KT');
        allEstAsrObservation_KT = [allEstAsrObservation_KT;
                                   asrObservations];
                                      
        asrObservations = runModelForward(p_KTAttn,A_KTAttn,B_KTAttn,C,D,v0,actualEEGObservations,T, 'KTAttn');
        allEstAsrObservation_KTAttn = [allEstAsrObservation_KTAttn;
                                       asrObservations];
    end
    
end

save(sprintf('latestResults_%d.mat', randi(1000)), 'allActualAsrObservation', 'allEstAsrObservation_KT', 'allEstAsrObservation_KTAttn');

%% Plot ROC curve

[X_KT,Y_KT] = perfcurve(allActualAsrObservation,allEstAsrObservation_KT,2);
[X_KTAttn,Y_KTAttn] = perfcurve(allActualAsrObservation,allEstAsrObservation_KTAttn,2);
figure; hold on;
    plot(X_KT,Y_KT, 'k--')
    plot(X_KTAttn,Y_KTAttn, 'k')
    legend('KT', 'KT-Attn');
    xlabel('False positive rate'); 
    ylabel('True positive rate')