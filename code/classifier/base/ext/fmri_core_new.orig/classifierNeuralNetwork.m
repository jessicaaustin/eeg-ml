%
% Train a neural network classifier using cross validation within
% the training set to decide when to stop. Uses the NetLab package,
% stored in the Netlab subdirectory
%
% In:
% - examples
% - labels
% - optionally, a classifierParameters cell array, containing:
%   - nHidden - how many hidden units (if 0, train a linear network)
%   - nIter   - how many iterations to train for
%               (if 0, use "leave 1 example out" cross-validation
%               within the training set)
%   - optimiz - optimization algorithm used in training
%
% Out:
% - a trained neural network in the NetLab format (the rest of the
% framework knows how to deal with it)
%
% History
% - 24 Nov 2004 - fpereira - created
%

function [net] = classifierNeuralNetwork( varargin )

DEBUG = 0;

l = length(varargin);
if l < 1; help classifierNeuralNetwork; return; end

examples = varargin{1}; [nExamples,nFeatures] = size(examples);
labels   = varargin{2}; sortedLabels = unique(labels); nLabels = length(sortedLabels);

classifierParameters = {};
if l > 2
  classifierParameters = varargin{3};
end

nIn     = nFeatures;
nHidden = nLabels;
nOut    = nLabels;
outfunc = 'softmax'; % good setting for classification problems
optimiz = 'scg'; % scaled conjugate gradient
iterParams = [];
nIter   = 0; % if 0 -> use CV within the trainining set to find # iterations

k = length(classifierParameters);
if k > 0
  nHidden = classifierParameters{1};
  if k > 1
    nIter = classifierParameters{2};
    if k > 2
      iterParams = classifierParameters{3};
      if k > 3
	optimiz = classifierParameters{4};
      end
    end
  end
end

%
% Train it
% 

% normalize the data and switch labels to 1-of-N encoding
examples = normalize(examples);
labelsN  = OneOfNencoding(labels);

% actual training
if nHidden
  % train a MLP
  net = mlp( nIn, nHidden, nOut, outfunc );
  % HACK to pass the network type with the network
  net.ourNetType = 'MLP';
else
  % if nHidden = 0, train a Generalised Linear Model
  net = glm( nIn, nOut, outfunc );
  net.ourNetType = 'GLM';
end
fprintf('training a %s neural network\n',net.ourNetType);

if nIter > 0
  % number of iterations has been specified    
else
  % determine nIter via cross validation within the training set
  nIter = computeOptimalNiter(examples,labels,net,nHidden,nIter,iterParams,DEBUG);
end

options     = zeros(1,18);
options(1)  = 1;                 % Print out error values    
options(14) = nIter;  

if 1
  % use same optimization algorithm for MLP and GLM
  net = netopt( net, options, examples, labelsN, optimiz);
else
  if nHidden
    net = netopt( net, options, examples, labelsN, optimiz);
  else
    net = glmtrain( net, options, examples, labelsN );
  end
end


%
% Cross validation within training set to determine how many
% iterations to train the network for
%

function [nIter] = computeOptimalNiter( varargin )

l = length(varargin);
examples = varargin{1};
labels   = varargin{2};
net      = varargin{3};
nHidden  = varargin{4};
nIter    = varargin{5};
iterParams = varargin{6};
if l > 6; DEBUG = varargin{7}; else; DEBUG = 0; end

sortedLabels = unique(labels); nLabels = length(sortedLabels);

for l=1:nLabels
  label   = sortedLabels(l);

  if nIter == 0
    % leave examples out

    places     = find(labels == label);
    nblocks(l) = length(places);
    for p=1:nblocks(l)
      orgSets{l}{p} = [places(p) places(p)];
    end

  elseif nIter == -1
    % leave 1 block out - HACK, DON'T USE

    %% Creates a data structure, orgSets{#labels}, where each entry is a
    %% cell array of intervals of indices into the array of
    %% examples. Each interval contains the beginning and end of
    %% the group of examples.
    
    % for the original examples
    % find the indices of the first image in each
    % block/presentation for this condition
    places  = find(labels == label);
    lastp   = find(diff(places)>1);
    breakp  = lastp + 1;
    breakp  = [1; breakp]; % block beginnings
    lastp   = [lastp; length(places)];
    nblocks(l) = length(breakp);

    % now compute the index intervals in the example array
    % corresponding to each block of examples
    orgSets{l} = cell(nblocks(l),1);
    for p=1:nblocks(l)
      orgSets{l}{p} = [places(breakp(p)) places(lastp(p))];
    end
  else
    fprintf('error: nIter=%d is not supported\n',nIter);pause;return
  end
    
end; % for over labels

%% check that all labels have the same number of trials
if sum(diff(nblocks))
  fprintf('error: # of trials should be the same for all labels\n');
  return;
else
  nTrialsPerLabel = nblocks(1);
end

if DEBUG
  fprintf('blocks for all\n');
  for l=1:nLabels
    label = sortedLabels(l);
    fprintf('\tlabel %d\t\n',label);
    for p=1:nblocks(l)
      fprintf('\t\t[%d,%d]\n',orgSets{l}{p});
      fprintf('\n');
    end
  end
  fprintf('\n');
  pause
end

%% Create a structure that says which examples are to be used in each fold.

    
if isempty(iterParams)
  nFolds = nTrialsPerLabel;
else
  nFolds = iterParams(1);
end

% Identify the image numbers that are going to be used
% as train and test in each fold
trainImagesPerFold    = cell(nFolds,1);
trainIntervalsPerFold = cell(nFolds,1);
testImagesPerFold     = cell(nFolds,1);
testIntervalsPerFold  = cell(nFolds,1);

for k=1:nFolds
  
  % find the test intervals
  testImagesPerFold{k} = [];
  testIntervalsPerFold{k} = {};
  
  for l=1:nLabels
    ii = orgSets{l}{k}; % image interval
    testImagesPerFold{k} = [testImagesPerFold{k},ii(1):ii(2)];
    testIntervalsPerFold{k}{l} = ii;
  end
  
  % find the train intervals
  trainImagesPerFold{k} = [];
  trainIntervalsPerFold{k} = {};idx = 1;
  
  for ak=1:k-1
    for l=1:nLabels
      ii = orgSets{l}{ak}; % image interval
      trainImagesPerFold{k} = [trainImagesPerFold{k},ii(1):ii(2)];
      trainIntervalsPerFold{k}{idx} = ii; idx = idx + 1;
    end
  end
  
  for ak=k+1:nFolds
    for l=1:nLabels
      ii = orgSets{l}{ak}; % image interval
      trainImagesPerFold{k} = [trainImagesPerFold{k},ii(1):ii(2)];
      trainIntervalsPerFold{k}{idx} = ii; idx = idx + 1;
    end
  end
end; % for over folds

% DEBUG: print out intervals
if DEBUG
  for k=1:nFolds
    fprintf('For fold %d\n',k);
    
    fprintf('\tTest %d\n',length(testImagesPerFold{k})); fprintf('\t\t');
    fprintf('%d ',testImagesPerFold{k});
    fprintf('\n');
    
    fprintf('\tTrain %d\n',length(trainImagesPerFold{k})); fprintf('\t\t');
    fprintf('%d ',trainImagesPerFold{k});
    fprintf('\n');
  end
  pause
end
   
%% Run the cross validation
%% - the examples used for train/test in each fold were defined in 4

% run k-fold cross validation

nIterPerFold = zeros(1,nFolds);
errorPerFold = zeros(1,nFolds);
size(labels)
for k=1:nFolds
  fprintf('Testing over fold %d\n',k);
  % a) select examples for test fold and train folds
  trainExamples   = examples(trainImagesPerFold{k},:);
  trainLabels     = labels(trainImagesPerFold{k},:);
  trainLabels1ofN = OneOfNencoding(trainLabels);
  testExamples    = examples(testImagesPerFold{k},:);
  testLabels      = labels(testImagesPerFold{k},:);
  testLabels1ofN  = OneOfNencoding(testLabels);

  % set up training
  nIterPerBurst   = 10;
  nBursts         = 20;
  errorAfterBurst{k} = zeros(1,nBursts);
  
  % train the network in little bursts
  foldNet = net;
  options = zeros(1,18);
  options(1)  = 1;                 % Print out error values
  options(14) = nIterPerBurst;
  method = 'scg';
  for b = 1:nBursts
    % train for <iterationsPerBurst>

    foldNet = netopt(foldNet, options, trainExamples, trainLabels1ofN, method);
    % apply to test set
    if nHidden
      % this is a MLP
      yt = mlpfwd(foldNet, testExamples);
    else
      % this is a GLM
      yt = glmfwd(foldNet, testExamples);
    end
      
    [yvalue,ypos]   = max(yt,[],2);
    predictedLabels = sortedLabels(ypos);
    errorAfterBurst{k}(b) = sum(predictedLabels~=testLabels)/length(testLabels);
  end

  [errorPerFold(k),nBurstsPerFold(k)] = min(errorAfterBurst{k});
end

nIter = ceil(median(nBurstsPerFold))*nIterPerBurst;

if DEBUG
  fprintf('error per fold:\t'); fprintf('%1.2f ',errorPerFold);fprintf('\n');
  fprintf('#bursts per fold:\t');fprintf('%d ',nBurstsPerFold);fprintf('\n');
  nIter
  pause
end
  
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

