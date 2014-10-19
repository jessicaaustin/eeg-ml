%
% Test a neural network classifier
%

function [] = classifierNeuralNetwork_test()

rand('seed',1685);
nExamplesPerClass = 20;
nFeatures         = 15000;
nClasses          = 3;
%nHidden           = nClasses; % train MLP
nHidden           = 0; % train GLM

[trainExamples,trainLabels] = generate3class(nExamplesPerClass,nFeatures);
[testExamples, testLabels]  = generate3class(nExamplesPerClass,nFeatures);
sortedLabels = [1 2 3];

%trainExamples = normalize(trainExamples);
testExamples  = normalize(testExamples);

if 0
  % test with direct calls to the routines
  
  % trained with # iterations selected by CV in the training set
  [net] = classifierNeuralNetwork(trainExamples,trainLabels,{nHidden});
  
  % trained with fixed # iterations
  %[net] = classifierNeuralNetwork(trainExamples,trainLabels,{nHidden,100});
  
  if nHidden
    yt = mlpfwd(net, testExamples);
  else
    yt = glmfwd(net, testExamples);
  end
  
  [yvalue,ypos]   = max(yt,[],2);
  predictedLabels = sortedLabels(ypos);
else
  % test using framework code
  models = trainClassifier(trainExamples,trainLabels,'neural',{nHidden});
  scores = applyClassifier(testExamples,models,'neural',{nHidden});
  [yvalue,ypos]   = max(scores,[],2);  
  predictedLabels = sortedLabels(ypos);
end


error = sum(predictedLabels~=testLabels')/length(testLabels)



%
% Ancillary code 
%

%% Transform a list of labels into a "1 of N" encoding
%% (a binary matrix with as many rows as examples, and as many
%% columns as labels. The row for an example is 1 in the column
%% corresponding to its label, and 0 everywhere else. The columns
%% are in the order of the labels sorted by value.

function [labels1ofN] = OneOfNencoding(labels)

classes   = unique(labels); nClasses = length(classes);
nExamples = length(labels);

labels1ofN = zeros(nExamples,nClasses);
for c = 1:nClasses
  label           = classes(c);
  labels1ofN(:,c) = (labels == label);
end


%% Normalize each feature to have mean 0 and standard deviation 1

function [Y] = normalize(X)

[nExamples,nFeatures] = size(X);
meanX = mean(X,1);
stdvX = std(X,0,1);

Y = X -  repmat(meanX,[nExamples,1]);
Y = Y ./ repmat(stdvX,[nExamples,1]);





function [examples,labels] = generate3class(nExamplesPerClass,nFeatures)

% Same, but now we have three classes
% 

% generate dataset of features with 0.5 bayes error

examples1 = zeros(nExamplesPerClass*3,nFeatures);

% 1/3 small error
sidx = 1; eidx = sidx + nFeatures/3 - 1;
examples1(1:nExamplesPerClass,sidx:eidx)                     = randn(nExamplesPerClass,nFeatures/3);
examples1(nExamplesPerClass+1:2*nExamplesPerClass,sidx:eidx) = randn(nExamplesPerClass,nFeatures/3) + 3;
examples1(2*nExamplesPerClass+1:end,sidx:eidx)               = randn(nExamplesPerClass,nFeatures/3) + 6;

% 1/3 medium error
sidx = eidx+1; eidx = sidx + nFeatures/3 - 1;
examples1(1:nExamplesPerClass,sidx:eidx)                     = randn(nExamplesPerClass,nFeatures/3);
examples1(nExamplesPerClass+1:2*nExamplesPerClass,sidx:eidx) = randn(nExamplesPerClass,nFeatures/3) + 1;
examples1(2*nExamplesPerClass+1:end,sidx:eidx)               = randn(nExamplesPerClass,nFeatures/3) + 2;

% 1/3 large error
sidx = eidx+1; eidx = sidx + nFeatures/3 - 1;
examples1(1:nExamplesPerClass,sidx:eidx)                     = randn(nExamplesPerClass,nFeatures/3);
examples1(nExamplesPerClass+1:2*nExamplesPerClass,sidx:eidx) = randn(nExamplesPerClass,nFeatures/3) + 0.05;
examples1(2*nExamplesPerClass+1:end,sidx:eidx)               = randn(nExamplesPerClass,nFeatures/3) + 0.1;

labels1   = ones(nExamplesPerClass*3,1)*1;
labels1(nExamplesPerClass+1:2*nExamplesPerClass,1) = 2;
labels1(2*nExamplesPerClass+1:3*nExamplesPerClass,1) = 3;

examples2 = zeros(nExamplesPerClass*3,nFeatures);

% 1/3 small error
sidx = 1; eidx = sidx + nFeatures/3 - 1;
examples2(1:nExamplesPerClass,sidx:eidx)                     = randn(nExamplesPerClass,nFeatures/3);
examples2(nExamplesPerClass+1:2*nExamplesPerClass,sidx:eidx) = randn(nExamplesPerClass,nFeatures/3) + 3;
examples2(2*nExamplesPerClass+1:end,sidx:eidx)               = randn(nExamplesPerClass,nFeatures/3) + 6;

% 1/3 medium error
sidx = eidx+1; eidx = sidx + nFeatures/3 - 1;
examples2(1:nExamplesPerClass,sidx:eidx)                     = randn(nExamplesPerClass,nFeatures/3);
examples2(nExamplesPerClass+1:2*nExamplesPerClass,sidx:eidx) = randn(nExamplesPerClass,nFeatures/3) + 1;
examples2(2*nExamplesPerClass+1:end,sidx:eidx)               = randn(nExamplesPerClass,nFeatures/3) + 2;

% 1/3 large error
sidx = eidx+1; eidx = sidx + nFeatures/3 - 1;
examples2(1:nExamplesPerClass,sidx:eidx)                     = randn(nExamplesPerClass,nFeatures/3);
examples2(nExamplesPerClass+1:2*nExamplesPerClass,sidx:eidx) = randn(nExamplesPerClass,nFeatures/3) + 0.05;
examples2(2*nExamplesPerClass+1:end,sidx:eidx)               = randn(nExamplesPerClass,nFeatures/3) + 0.1;

labels2   = ones(nExamplesPerClass*3,1)*1;
labels2(nExamplesPerClass+1:2*nExamplesPerClass,1) = 2;
labels2(2*nExamplesPerClass+1:3*nExamplesPerClass,1) = 3;

examples = [examples1;examples2];
labels   = [labels1; labels2];
