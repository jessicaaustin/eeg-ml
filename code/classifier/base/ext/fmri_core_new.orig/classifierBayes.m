% Train several flavours of bayesian classifier
%
% [models] = classifierBayes( examples, labels, [variation],[parameters] )
%
% In:
% - training Examples and Labels
% - optionally
%   - variation - defaults to 'nbayes', a vanilla GNB, but can also be
%     - nbayesPooled - GNB with variance estimates pooled across classes
%     - lda          - gaussian classifier with pooled covariance matrix
%     - ldaSVD       - LDA after SVD dimensionality reduction
%
% Out:
% - a models structure containing the classifier
%   (please see comments for trainClassifier_v2 for more details)
%
% Dependencies
%
% History:
% - 28 Oct 04 - fp - created from classifierNBayes code
% - 3  Jun 05 - fp - fixed a bug in reading classifierParameters
%
% Notes:
%

function [models] = classifierBayes( varargin )

%
% process arguments
%

l = length(varargin);
if l < 2; help classifierBayes; return; else

  if l > 2
    classifier = varargin{3};
    classifierParameters = {};
    if l > 3; classifierParameters = varargin{4}; end
  else
    % vanilla naive bayes
    classifier = 'nbayes';
  end 

  switch classifier
   case {'nbayes','nbayesPooled','nbayes-unitvariance','lda','qda','ldaSVD','ldaCV','qdaSVD','qdaCV','nbayesPooledResampling'}
    % we know about these
   otherwise
    fprintf('classifierNBayes: I don''t know about classifier %s\n',classifier);
  end

end

trainSet    = varargin{1};
trainLabels = varargin{2};

% some methods require processing of the data that affects the
% number of dimensions, so do it before all that is figured out

switch classifier
 case {'ldaSVD','ldaCV','qdaSVD','qdaCV'}
  % the training set returned has as many dimensions as components
  % in the basis of eigenimages.
  
  % Each row in trainSet is now the coordinates of that example on
  % the basis of eigenimages that explains <fractionToKeep> of the
  % variance
  fractionToKeep = 0.9;
  if ~isempty(classifierParameters); fractionToKeep = classifierParameters{1};end

  switch classifier
   case {'ldaSVD','qdaSVD'}
    [trainSet,V] = reduceDimensionality(trainSet,trainLabels,'SVD',fractionToKeep);
   case {'ldaCV','qdaCV'}
    % first use SVD to reduce dimensionality (canvar code is too inefficient)
    [trainSet,Vsvd] = reduceDimensionality(trainSet,trainLabels,'SVD',fractionToKeep);
    [trainSet,Vcv]  = reduceDimensionality(trainSet,trainLabels,'CV',fractionToKeep);
  end  
  
  % crop?
  
 otherwise
end

nTrain      = size(trainSet,1);
nFeatures   = size(trainSet,2);
nLabels     = size(trainLabels,2);
  
sortedLabelValues = sort(unique(trainLabels));
nClasses          = length(sortedLabelValues);

weights   = zeros(nFeatures+1,nClasses);

models = cell(nClasses+3,1); % will keep the final results

%
% Train
%

% training involves learning means for each feature and either
% standard deviations for each feature or a full covariance matrix

% these are stored in these two matrices, each of dim (#classes,#features)
%% if space is not preallocated each new element stored forces a reallocation
means = zeros(nClasses,nFeatures);
stds  = zeros(nClasses,nFeatures);
  
% training also involves learning classPriors (#classes)
classPriors = zeros(nClasses,1);
nPerClass   = zeros(nClasses,1);


%% a) Find the examples belonging to each class and their means/stdevs

indices = cell(nClasses,1);

% hack - save it now, as the training data will have its mean subtracted
switch classifier
 case {'nbayesPooledResampling'}
  models{nClasses+3}{1} = trainSet;
  models{nClasses+3}{2} = trainLabels;
end


for c = 1:nClasses
  cond = sortedLabelValues(c);
  indices{c}   = find( trainLabels == cond );
  nPerClass(c) = length(indices{c});

  % compute means and stdevs per class
  means(c,:) = mean(trainSet(indices{c},:),1);

  switch classifier
   case {'nbayes'}
    stds(c,:)  = std(trainSet(indices{c},:),0,1);
    
   case {'nbayes-unitvariance'}
    stds(c,:)  = ones(1,nFeatures);
    
   case {'nbayesPooled','nbayesPooledResampling','lda','qda','ldaSVD','ldaCV','qdaSVD','qdaCV'}
    % subtract means from data, as it will help
    % compute pooled feature standard deviations
    % or a full covariance matrix
    
    trainSet(indices{c},:) = trainSet(indices{c},:) - repmat(means(c,:),nPerClass(c),1);
  end
end

classPriors = nPerClass / nTrain;

%% b) Standard deviations/Cov matrix for other classifiers

switch classifier
 case {'nbayes'}
  % all done
  
 case {'nbayesPooled','nbayesPooledResampling'}
  % estimate standard deviations over the entire data, now that
  % they have a common mean of 0 in each feature
  stds = repmat( std(trainSet,0,1), nClasses,1 );

 case {'lda','ldaSVD','ldaCV'}
  % estimate a full covariance matrix for the centred training set
  if isempty(classifierParameters) useRobustEstimate = 0; else
  useRobustEstimate = classifierParameters{1}; end
  fprintf('apply_classifier: useRobustEstimate = %d\n',useRobustEstimate);
  
  if ~useRobustEstimate
    stds = cov(trainSet);
  else
    % use estimate from LIBRA
    rew  = mcdcov(trainSet,'plots',0);
    stds = rew.cov; 
  end    
    
 case {'qda','qdaSVD','qdaCV'}
  % estimate one covariance matrix per class
  if isempty(classifierParameters) useRobustEstimate = 0; else
  useRobustEstimate = classifierParameters{1}; end
  fprintf('apply_classifier: useRobustEstimate = %d\n',useRobustEstimate);
  
  classCovarianceMatrices = cell(nClasses,1);
  for c = 1:nClasses

    if ~useRobustEstimate      
      classCovarianceMatrices{c} = cov(trainSet(indices{c},:));
    else
      rew = mcdcov(trainSet(indices{c},:),'plots',0);
      classCovarianceMatrices{c} = rew.cov; 
    end
  end
end
  

%% c) pack the results into a cell array



% Generative model - each cell contains a cell array
% where each cell has one parameter - mean, covariance matrix, etc
for c=1:1:nClasses
  models{c} = cell(2,1);
  models{c}{1} = means(c,:);
  switch classifier
   case {'nbayes','nbayesPooled','nbayes-unitvariance','nbayesPooledResampling'}
    models{c}{2} = stds(c,:);    
   otherwise
    models{c}{2} = [];
  end
  models{c}{3} = classPriors(c);
end

% Discriminative model - a cell array of sets of weights
models{nClasses+1} = cell(1,1);
  
% Training Set information
trainingSetInfo.nExamples         = nTrain;
trainingSetInfo.nFeatures         = nFeatures;
trainingSetInfo.nClasses          = nClasses; % same thing as labels
trainingSetInfo.sortedLabelValues = sortedLabelValues;
trainingSetInfo.classPriors       = classPriors;
models{nClasses+2} = trainingSetInfo;

% Extra information (depends on classifier)

switch classifier
 case {'nbayes','nbayesPooled','nbayes-unitvariance'}
  % nothing
 case {'nbayesPooledResampling'}
  % the training data has been put here already
 case {'lda'}
  models{nClasses+3}{1} = stds; % covariance matrix
 case {'qda'}
  models{nClasses+3}{1} = classCovarianceMatrices;
 case {'ldaSVD','ldaCV'}
  models{nClasses+3}{1} = stds; % covariance matrix
  if isequal(classifier,'ldaSVD')
    models{nClasses+3}{2} = V; % projection into SVD space
  else
    models{nClasses+3}{2}{1} = Vsvd; % project into SVD space
    models{nClasses+3}{2}{2} = Vcv;  % project into CV space of that
  end
  
 case {'qdaSVD','qdaCV'}
  models{nClasses+3}{1} = classCovarianceMatrices; % covariance matrix
  if isequal(classifier,'qdaSVD')
    models{nClasses+3}{2} = V; % projection into SVD space
  else
    models{nClasses+3}{2}{1} = Vsvd; % project into SVD space
    models{nClasses+3}{2}{2} = Vcv;  % project into CV space of that
  end
end
	
% Global model - trained over the entire dataset
%models{nClasses+3} = cell(2,1);
%models{nClasses+3}{1} = mean(trainSet,1);
%models{nClasses+3}{2} = std(trainSet,0,1);



%% reduce dimensionality of the training set

function [MM,V] = reduceDimensionality(trainSet,trainLabels,method,fraction)

switch method

  case {'SVD'}
   % SVD dataset
   [U,S,V] = compute_fastSVD(trainSet);
   [nExamples,nVoxels] = size(trainSet);  
   maxp = min(nExamples,nVoxels);
   r = rank(S);
   
   % reduce the three matrices
   U = U(:,1:maxp);
   if nExamples < nVoxels; S = S(:,1:maxp); else; S = S(1:maxp,:); end 
   V = V(:,1:maxp);
   
   % scree plot intrinsic dimensionality
   scores           = diag(S);
   normalizedScores = scores.^2;
   normalizedScores = normalizedScores ./ sum(normalizedScores);
   components       = V';
   MM         = U*S;
   
   % find those components - sort and invert order (larger first)
   [sns,sindices] = sort(normalizedScores);
   sns      = flipud(sns);
   sindices = flipud(sindices);
   
   % Swap components into the appropriate order
   components = components(sindices,:);
   scores     = scores(sindices);
   normalizedScores = normalizedScores(sindices);  
   MM         = MM(:,sindices);

   % find the components to keep in order to account for <fraction> of
   % the variance
   tmp = cumsum(sns);
   pos = find(tmp >= fraction);
   nComponentsToKeep = pos(1); % minimum number of components;
   nComponentsToKeep = min(nComponentsToKeep,r);
   
   sindices = 1:nComponentsToKeep;
   
   % crop the matrices accordingly
   components = components(sindices,:);
   scores     = scores(sindices);
   normalizedScores = normalizedScores(sindices);  
   MM         = MM(:,sindices);
   V = components';  
   
   nComponents = length(sindices);
  
 case {'CV'}
  
  trainLabelsOneOfN = OneOfNencoding(trainLabels);  
%  [cvals,V] = compute_canvar( trainSet, trainLabelsOneOfN );
  [cvals,V] = canvar( trainSet, trainLabelsOneOfN );
  MM = trainSet * V;
  
end


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



function test_2D()

%% Shared by all tests
% nPerclass*2 examples, nPerClass for class A, then nPerClass for B,
% train will be first nPerClass/2 of each class, test the remainder

nPerClass = 100;
nTotal    = 2*nPerClass;
trainRange = [1:nTotal/4,nTotal/2+1:3*nTotal/4];
testRange  = [nTotal/4+1:nTotal/2,3*nTotal/4+1:nTotal];
labels     = [repmat(1,nPerClass,1);repmat(2,nPerClass,1)];
labelsTrain   = labels(trainRange);
labelsTest    = labels(testRange);

%% NB assumption holds
sigma = [[1 0];[0 1]];

% means are apart
mu1   = [0 0];
mu2   = [2 0];
a = mvnrnd(mu1,sigma,nPerClass);
b = mvnrnd(mu2,sigma,nPerClass);
examples   = [a;b];
examplesTrain = examples(trainRange,:);
examplesTest  = examples(testRange,:);

% means are close 
mu1   = [0 0];
mu2   = [0.1 0];
a = mvnrnd(mu1,sigma,nPerClass);
b = mvnrnd(mu2,sigma,nPerClass);
examples   = [a;b];
examplesTrain = examples(trainRange,:);
examplesTest  = examples(testRange,:);


%% NB assumption doesn't hold
sigma = [[1 1];[1 1]];

% means are apart
mu1   = [0 0];
mu2   = [2 0];
a = mvnrnd(mu1,sigma,nPerClass);
b = mvnrnd(mu2,sigma,nPerClass);
examples   = [a;b];
examplesTrain = examples(trainRange,:);
examplesTest  = examples(testRange,:);

classifier = 'nbayesPooledResampling';
[models] = classifierNBayes_v2(examplesTrain,labelsTrain,classifier);
[scores] = applyClassifier_v2(examplesTest,models,classifier);
[result1,learntLabels] = summarizePredictions_v2(examplesTest,models,scores,classifier,'accuracy',labelsTest);

classifier = 'nbayesPooled';
[models] = classifierNBayes_v2(examplesTrain,labelsTrain,classifier);
[scores] = applyClassifier_v2(examplesTest,models,classifier);
[result2,learntLabels] = summarizePredictions_v2(examplesTest,models,scores,classifier,'accuracy',labelsTest);


[models] = classifierNBayes_v2(examplesTrain,labelsTrain,'lda');
[scores] = applyClassifier_v2(examplesTest,models,'lda');
[result,learntLabels] = summarizePredictions_v2(examplesTest,models,scores,'nbayes','accuracy',labelsTest);

[models] = classifierNBayes_v2(examplesTrain,labelsTrain,'ldaSVD');
[scores] = applyClassifier_v2(examplesTest,models,'ldaSVD');
[result,learntLabels] = summarizePredictions_v2(examplesTest,models,scores,'nbayes','accuracy',labelsTest);


models{1}{1}
models{1}{2}
models{2}{1}
models{2}{2}




clf;
hold on;
plot(a(:,1),a(:,2),'b.','MarkerSize',6);
plot(b(:,1),b(:,2),'r.','MarkerSize',6);
hold off;



%
%


  
function test_simple()

  % simple, balanced normal dataset (!= means)
  
  % each class generated from a 10-D normal
  means{1} = 1:1:10;
  means{2} = 2 + means{1};
  examples = zeros(40,10);
  labels   = zeros(40,1);
  
  for i=1:1:40
    if mod(i,2)
      examples(i,:) = randn(1,10) + means{1};
      labels(i) = 2;
    else
      examples(i,:) = randn(1,10) + means{2};
      labels(i) = 1;
    end
  end
  
  expInfo.experiment = 'simple';
  expInfo.meta       = [];
  
  trainSet    = examples(1:30,:);
  testSet     = examples(31:40,:);
  trainLabels = labels(1:30);
  testLabels  = labels(31:40);
  
  [models] = classifierLDA(trainSet, trainLabels);
  [scores] = applyClassifier_v2(testSet,models,'nbayes');
  [result,learntLabels] = summarizePredictions(testSet,models,scores,'nbayes','accuracy',testLabels);
  result{1}
  
  % check models{2} {1} and {2} for the mean and variance of class 0
  
function test_unbalanced()

  % unbalanced: 1 of 0 to 5 of 1
  
  % each class generated from a 10-D normal
  means{1} = 1:1:10;
  means{2} = 2 + means{1};
  examples = zeros(84,10);
  labels   = zeros(84,1)
  
  % train set
  for i=1:1:50
    examples(i,:) = randn(1,10) + means{1};
    labels(i) = 2;
  end
  
  for i=51:1:60
    examples(i,:) = randn(1,10) + means{2};
    labels(i) = 1;
  end
    
  % test set
  for i=61:1:80
    examples(i,:) = randn(1,10) + means{1};
    labels(i) = 2;
  end
  
  for i=81:1:84
    examples(i,:) = randn(1,10) + means{2};
    labels(i) = 1;
  end

  % pick set
  trainSet    = examples(1:60,:);
  testSet     = examples(61:84,:);
  trainLabels = labels(1:60,:);
  testLabels  = labels(61:84,:);

  [models] = classifierNBayes(trainSet, trainLabels);
  [scores] = applyClassifier(testSet,models,'nbayes');
  [result,learntLabels] = summarizePredictions(testSet,models,scores,'nbayes','accuracy',testLabels);
  result{1}

  % check models{2} {1} and {2} for the mean and variance of class
  % 0

  
function test_multiple()

  % each class generated from a 10-D normal
  means{1} = 1:1:10;
  means{2} = 2 + means{1};
  means{3} = 4 + means{1};
  means{4} = 6 + means{1};
  means{5} = 8 + means{1};
  means{6} = 10 + means{1};
  examples = zeros(6*(10+2),10);
  labels = zeros(6*(10+2),1);

  expInfo.experiment = 'simple';
  expInfo.meta       = [];
  
  % train set
  idx = 1;
  for i=1:1:10
    for c=1:1:6
      examples(idx,:) = randn(1,10) + means{c};
      labels(idx) = c;
      idx = idx + 1;
    end
  end
  nTrain = 10*6;
  
  % test set
  for i=1:1:2
    for c=1:1:6
      examples(idx,:) = randn(1,10) + means{c};
      labels(idx) = c;
      idx = idx + 1;
    end
  end
  nTest = 2*6;
  
  % pick set
  trainSet    = examples(1:nTrain,:);
  testSet     = examples((nTrain+1):(nTrain+nTest),:);
  trainLabels = labels(1:nTrain,:);
  testLabels  = labels((nTrain+1):(nTrain+nTest),:);
  
  [models] = classifierNBayes(trainSet, trainLabels);
  [scores] = applyClassifier(testSet,models,'nbayes');  
  [result,learntLabels] = summarizePredictions(testSet,models,scores,'nbayes','accuracy',testLabels);
  result{1}
  % check models{2} {1} and {2} for the mean and variance of class
  % 0

  
  [models,scores] = trainClassifierL1O(examples,labels,'nbayes','full');
  [result,learntLabels] = summarizePredictions(examples,models,scores,'nbayes','accuracy',labels);

  [models1,scores1] = trainClassifier_kFoldCV(examples,labels,expInfo,'nbayes','full',72);
  [result1,learntLabels1] = summarizePredictions(examples,models1,scores1,'nbayes','accuracy',labels);
