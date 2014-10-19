%
% A simple test
%
% For each classifier, train and test on a sample dataset


function [] = run_testSuite1()

%
% Classifiers to use
% 

DEBUG = 0;

if ~DEBUG
  classifiers = {'nbayes','nbayesPooled','nbayes-unitvariance','svm','knn','neural','neural','logisticRegression','pairwise'};

  c = 1;
  classifierParameters{c} = {};    c=c+1;
  classifierParameters{c} = {};    c=c+1;
  classifierParameters{c} = {};    c=c+1;
  classifierParameters{c} = {};    c=c+1;
  classifierParameters{c} = {};    c=c+1;
  classifierParameters{c} = {1,1000}; c=c+1; % 1 hidden unit, at most 1000 iterations
  classifierParameters{c} = {2,1000}; c=c+1; % 2 hidden units,leave-1-out stopping
  classifierParameters{c} = {};    c=c+1;
  classifierParameters{c} = {0,0}; c=c+1; % linear units only,l-1-out stopping
else
  % old classifiers
%  classifiers = {'svm','nnets','nbayes-unitvariance'}
  classifiers = {'nbayes-unitvariance'}
c = 1;
%  classifierParameters{c} = {};    c=c+1;
%  classifierParameters{c} = {};    c=c+1;
  classifierParameters{c} = {};    c=c+1;  
end

%
% Create a dataset
%

nFeaturesTotal = 1000; % how many features total
nFeaturesGood  = 20; % how many of those are not noise
meanDifference = 1; % how far apart are class means in those
featureStdev   = 1; % what is the standard deviation of the class densities
nPerClass      = 20; % how many examples to generate per class


trainingSet = randn(nPerClass*2,nFeaturesTotal) * featureStdev;
trainingSet(1:nPerClass,1:nFeaturesGood) = trainingSet(1:nPerClass,1:nFeaturesGood) + meanDifference;

trainingLabels = repmat(2,nPerClass*2,1); trainingLabels(1:nPerClass) = 1;
testingLabels  = trainingLabels;

testingSet = randn(nPerClass*2,nFeaturesTotal) * featureStdev;
testingSet(1:nPerClass,1:nFeaturesGood) = testingSet(1:nPerClass,1:nFeaturesGood) + meanDifference;

%
% Run the test
%

ofd = fopen('testSuite1.results','w');
classifers_accuracy=[];

for c = 1:7%length(classifiers)

  classifier       = classifiers{c};
  classifierParams = classifierParameters{c}; 
  
  [trainedClassifier] = trainClassifier(trainingSet,trainingLabels,classifier,classifierParams);
  trainedClassifier
  [scores] = applyClassifier(testingSet,trainedClassifier);
  [result] = summarizePredictions(scores,trainedClassifier,'accuracy',testingLabels);
  
  classifers_accuracy = [classifers_accuracy result{1}];
  fprintf(ofd,'%s\t%1.2f\n',classifier,result{1});
  fprintf('%s\t%1.2f\n',classifier,result{1});
end
classifers_accuracy

%
% Match classifier name to corresponding function
%

function [classifierFunction] = classifierNameToFunction( classifier )

switch classifier
  
 case {'nbayes'}
  
 case {'nbayesPooled'}
  
 case {'knn'}
  
 case {'svm'}
  
 case {'neuralNetwork'}
  
 case {'qda'}
  
 case {'lda'}

end
