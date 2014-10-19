% working example illustrating use of fmri_core software to load
% starplus data, transform it into examples, train classifier,
% apply it to test data, and see the result.
%
% Note to execute this you must be on a machine containing the data
% for the 'data-syntamb2' study for subject 02945_22.
%
% - April 13, 2005 created by Wei Wang
%
% - Test: April 13, 2005
%

subject = '02945_22';
rois    = {'LB' 'LT'};
study = 'data-syntamb2';

% do THIS
% returns cell arrays of IDMs, one per ROI
separateRois = 1;
[is,ds,ms] = loadSubjectdataMult(study,subject,rois,separateRois);
% use LB
info = is{1}; data = ds{1}; meta = ms{1};

% or THIS - gives LB+LT merged together
[info,data,meta] = loadSubjectdata(study,subject,rois);

% pick active voxels
[ainfo,adata,ameta] = transformIDM_selectActiveVoxels(info,data,meta,5);

% average all trials in each condition
[avgainfo,avgadata,avgameta] = transformIDM_avgTrialCondition(ainfo,adata,ameta);

% plot voxel in column 3
plotVoxel(avgainfo,avgadata,avgameta,3);

% transform one of them into examples
[examples,labels,expInfo] = idmToExamples_fixation( is{1},ds{1},ms{1});

% split the data in half
ntotal   = size(examples,1);
oindices = 1:2:(ntotal-1);
eindices = 2:2:ntotal;
trainExamples = examples(oindices,:);
trainLabels   = labels(oindices,1);
testExamples  = examples(eindices,:);
testLabels    = labels(eindices,1);

% train a classifier
classifiers = {'nbayes','nbayesPooled','nbayes-unitvariance','svm','knn','neural','neural','logisticRegression','SMLR'};

c = 1;
classifierParameters{c} = {};    c=c+1;
classifierParameters{c} = {};    c=c+1;
classifierParameters{c} = {};    c=c+1;
classifierParameters{c} = {};    c=c+1;
classifierParameters{c} = {};    c=c+1;
classifierParameters{c} = {1,1000}; c=c+1; % linear network, at most 1000 iterations
classifierParameters{c} = {2,1000}; c=c+1; % 2 hidden units,at most 1000 iterations
classifierParameters{c} = {};    c=c+1;
classifierParameters{c} = {}; c=c+1; 

ofd = fopen('example_syntamb2.results','w');
classifers_accuracy=[];

for c = 1:length(classifiers)

  classifier       = classifiers{c};
  classifierParams = classifierParameters{c}; 
  
  [trainedClassifier] = trainClassifier(trainExamples,trainLabels,classifier,classifierParams);
  trainedClassifier
  [scores] = applyClassifier(testExamples,trainedClassifier);
  [result] = summarizePredictions(scores,trainedClassifier,'accuracy',testLabels);
  
  classifers_accuracy = [classifers_accuracy result{1}];
  fprintf(ofd,'%s\t%1.2f\n',classifier,result{1});
  fprintf('%s\t%1.2f\n',classifier,result{1});
end

classifers_accuracy
