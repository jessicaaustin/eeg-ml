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
        actualEEGObservations = thresholdAndFillAttention(sequences.attention{si}, sequences.timeelapsed{si}, sequences);
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
hFig=figure; hold on;
    plot(X_KT,Y_KT, 'b', 'LineWidth', 3)
    plot(X_KTAttn,Y_KTAttn, 'r--', 'LineWidth', 3)
    l=legend('KT', 'KT-Attn', 'Location', 'southeast');
    xlabel('False positive rate'); 
    ylabel('True positive rate')
    
    figureHandle = gcf;
    set(findall(figureHandle,'type','text'),'fontSize',14);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',14)
        
    set(hFig, 'Position', [0 0  400 300])

    
%% Plot results

load('latestResults_KTAttnSlightlyBetter_NoTimeElapsed');
allEstAsrObservation_KTAttn_NoTE = allEstAsrObservation_KTAttn;

load('latestResults_KTAttnSlightlyBetter_ThresholdTimeElapsedPos');
allEstAsrObservation_KTAttn_TE = allEstAsrObservation_KTAttn;

[X_KT,Y_KT] = perfcurve(allActualAsrObservation,allEstAsrObservation_KT_3,2);
[X_KTAttn_NoTE,Y_KTAttn_NoTE] = perfcurve(allActualAsrObservation,allEstAsrObservation_KTAttn_NoTE,2);
[X_KTAttn_TE,Y_KTAttn_TE] = perfcurve(allActualAsrObservation,allEstAsrObservation_KTAttn_TE,2);
hFig=figure; hold on;
    plot(X_KT,Y_KT, 'b', 'LineWidth', 3)
    plot(X_KTAttn_NoTE,Y_KTAttn_NoTE, 'm--', 'LineWidth', 3)
    plot(X_KTAttn_TE,Y_KTAttn_TE, 'r--', 'LineWidth', 3)
    l=legend('KT', 'KT-Attn (EEG only)', 'KT-Attn (EEG with time elapsed)', 'Location', 'southeast');
    xlabel('False positive rate'); 
    ylabel('True positive rate')
    
    figureHandle = gcf;
    set(findall(figureHandle,'type','text'),'fontSize',14);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',14)
        
    set(hFig, 'Position', [0 0  600 500])

    