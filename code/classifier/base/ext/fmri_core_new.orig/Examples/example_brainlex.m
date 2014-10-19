% working example illustrating use of fmri_core software to visualize, load
% starplus data, transform it into examples, train classifier,
% apply it to test data, and see the result.
%
% Note to execute this you must be on a machine containing the data
% for the 'data-brainlex' study for subject 08057.
%
% - April 13, 2005 created by Wei Wang
%
% - Test: April 13, 2005
%


study   = 'data-brainlex';
subject = '08057';
rois={'CALC' 'LT'};

% returns cell arrays of IDMs
[i,d,m] = loadSubjectdata(study,subject,rois);

% plot an image showing the brain activity for the 7th snapshot in trial 3,
% displaying the z=10 slice
n=7; trialNum=3; z=10; 
plotSnapshot(i,d,m,trialNum,z,n,0,0);

% watch a movie of the brain activity for one z slice of trial 3
z=10
M=animateTrial(i,d,m,3,z);

% watch it again
movie(M);

% watch a movie of the brain activity for all 16 z slices of
% trial 3
M=animate16Trial(i,d,m,3);


% erase the big data structure
clear M;

% transform one IDM them into examples
[examples,labels,expInfo] = idmToExamples_fixation(i,d,m,'full');

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

ofd = fopen('example_brainlex.results','w');
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
