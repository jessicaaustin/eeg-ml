% working example illustrating use of fmri_core software to load
% starplus data, transform it into examples, train classifier,
% apply it to test data, and see the result.
%
% Note to execute this you must be on a machine containing the data
% for the 'data-starplus-sp' study for subject 04847_20_sp3terC.
%
% - April 13, 2005 created by Wei Wang
% - adapted from FP
%
% - Test: April 13, 2005
%


subject = '04847_20_sp3terC' 
rois    = {'CALC'};
study   = 'data-starplus-sp';

% load the data - returns IDM for one ROI
[info,data,meta] = loadSubjectdata(study,subject,rois);

% transform one IDM them into examples
[examples,labels,expInfo] = idmToExamples_fixation(info,data,meta,'full');
% examples is now a mxn array of examples; labels is now a nx1 column
% vector of labels,  where m is the number of examples, and n is the
% number of features

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
classifierParameters{c} = {1,'correlation'};    c=c+1;
classifierParameters{c} = {1,1000}; c=c+1; % linear network, at most 1000 iterations
classifierParameters{c} = {2,1000}; c=c+1; % 2 hidden units,at most 1000 iterations
classifierParameters{c} = {};    c=c+1;
classifierParameters{c} = {};    c=c+1;

ofd = fopen('example_starplus.results','w');
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

