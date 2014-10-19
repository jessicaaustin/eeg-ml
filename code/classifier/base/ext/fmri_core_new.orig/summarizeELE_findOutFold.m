% Rank examples by some measure of how much of an outlier they are.
% The methods in this function share the fact that they work using
% balanced leave 1 out,learning some model on the training set for
% each fold and using it to score the examples on the respective
% test set.
%
% In:
% - examples
% - labels
% - method:
%   - SVD - runs SVD on the training set, uses learnt basis to
%           reconstruct the examples in the test set.
%   - nbayesPooled|logisticRegression2
%         - trains one of these on the training set, applies the
%         classifier to the test set to generate a set of <#
%         classes> scores for each example, uses those scores to decide.
%   - textFile
%         <rank according to a given file - not implemented yet>
%   - random
%         - rank at random
%
% - method parameters, a cell array
%   - SVD parameters - {<subtract mean of data matrix>,<percentage of variance to keep>}
%   - nbayesPooled - {}
%   - logistic regression - {<lambda>,<step size>}
%
% - optionally (doesn't have to be specified)
%   - FS method to select features on each training set
%   - FS params
%   - # voxels to select
%
% - optionally
%   - crossValidationMethod (default is 'leave1out', leaving out one example of each class)
%
% Out:
% - a ranking of the examples, from most outlying to least outlying
%   (at the end of the function is code that knows whether the
%   score produced by each method is good or bad and ranks accordingly)
%
% Notes:
%
% - There's a toy example to try this with at the bottom of the code.
% - This code needs to know where to find at least  "summarizeELE_rankByCCBImethod.m"
%
% Examples:
%
% er=summarizeELE_findOutFold(examples,labels,'SVD',{0,0.9});
% er=summarizeELE_findOutFold(examples,labels,'SVD',{0,0.9},'CCBI',{'replicability'},100);
% er=summarizeELE_findOutFold(examples,labels,'nbayesPooled',{});
% er=summarizeELE_findOutFold(examples,labels,'logisticRegression2',{1 0.01},'CCBI',{'replicability'},100);
%
% History:
% - 2005 Sep 26 - fpereira - created

function [exampleRanking,exampleScores] = summarizeELE_findOutFold(varargin)

%
% Process parameters and figure things out
%
rand('state',1685);

DEBUG = 0;
this = 'summarizeELE_findOutFold';
l = length(varargin);
if l < 4
  eval(sprintf('help %s;',this));pause;return;
else
  examples = varargin{1};
  labels   = varargin{2};
  method   = varargin{3};
  params   = varargin{4};
  useFS = 0;
  nVoxelsToUse = 0;
  if l > 4
    if l == 7
      FSmethod     = varargin{5};
      FSparams     = varargin{6};
      nVoxelsToUse = varargin{7};
      useFS        = 1;
    else
      fprintf('%s: to use FS you must specify all parameters',this);pause;return;
    end
  end
  
  [nExamples,nFeatures] = size(examples);
  crossValidationMethod = 'leave1out';
  useClass = 0;
  if l > 7;   crossValidationMethod = varargin{8}; end

  % deal with method parameters
  switch method

   case {'SVD'}
    % parameters
    if length(params) < 2; fprintf('%s: SVD parameter error',this);pause;return; end
    percentageOfVariance = params{1};
    subtractFeatureMean  = params{2};
    
   case {'nbayesPooled','nbayes'}
    % no parameters required
  
   case {'logisticRegression2'}
    lambda   = 1;
    stepSize = 0.001;
    if length(params) > 0;   lambda   = params{1};
      if length(params) > 1; stepSize = params{2};
      end
    end

   case {'textFile'}
    % load the ranking from a text file
    textFile = params{1};
    load(textFile);
    % higher is better, so flip the ranking
    [exampleScores,exampleRanking] = sort(subjectScores);
    exampleRanking = flipud(exampleRanking);
    exampleScores  = flipud(exampleScores);
    
   case {'random'}
    % nothing to do
    
   otherwise
    eval(sprintf('help %s;',this));pause;return;  
  end
  
end

classes = unique(labels); nClasses = length(classes);
[nExamples,nFeatures] = size(examples);

if 0
  %% find the examples belonging to each class and compute statistics
  for c = 1:nClasses
    class  = classes(c);  
    examplesClass{c}    = find(labels == class);
    examplesInClass{c}  = examples(examplesClass{c},:);
    nExamplesClass(c)   = length(examplesClass{c});  
    meanExampleClass{c} = mean(examplesInClass{c},1);
  end  
  clear examples
end
  
%% example information about the processing steps performed on
%% the set of examples we start with, with respect to the
%% original set of examples
%% - examples
%% - features
examplesKeptMask = ones(size(examples,1),1);
originalLabels   = labels;
exampleInfo.examplesKeptMask    = examplesKeptMask;
exampleInfo.examplesKeptIndices = find(examplesKeptMask);
exampleInfo.originalLabels      = originalLabels;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Run classification loop (train + apply)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% figure out cross-validation

% Create the division of examples into folds,  by determining
% what indices(rows) in the example matrix match the blocks of
% examples for each condition

% how many successive examples of the same class go into a fold
% hardwired at defaults, for now
switch crossValidationMethod
 case {'leave1out','evenOdd','leave1outExample'}
  exampleBlockSize = 1;	 
 case {'leave1outUnbalanced'}
  exampleBlockSize = 0;
end

[foldExamplesTrain,foldExamplesTest] = aux_divideIntoFolds(labels,crossValidationMethod,exampleBlockSize);
nFolds = length(foldExamplesTrain);
  
%% run the classification loop

% assemble one array of scores per <# features> to try
scoresPerFold = cell(nFolds,1);
nTestPerFold  = zeros(1,nFolds);

for k = 1:nFolds

  fprintf('%s: processing fold %d of %d\n',this,k,nFolds);
      
  %% create train and test sets for this fold
  trainExamples = examples(foldExamplesTrain{k},:);
  trainLabels   = labels(foldExamplesTrain{k},:);
  testExamples  = examples(foldExamplesTest{k},:);
  testLabels    = labels(foldExamplesTest{k},:);

  if useFS
    % foldExamplesTrain{k} are indices with respect to examples
    % at the beginning of the fold loop, i.e. after outlier removal.
    % Find the indices taken with respect to the original examples.
    
    if 0
      exampleInfo.examplesKeptMask    = examplesKeptMask;
      exampleInfo.examplesKeptIndices = find(examplesKeptMask);
      exampleInfo.labelsOrg           = labels;
    end
    
    tmp = (1:length(examplesKeptMask))';
    tmp = tmp(exampleInfo.examplesKeptIndices);
    
    trainExamplesOrgIndices = tmp(foldExamplesTrain{k});
    testExamplesOrgIndices  = tmp(foldExamplesTest{k});
    
    %% rank features using the selected method
    [sortedFeatures,sortedValues] = aux_rankFeatures(trainExamples,trainLabels,FSmethod,FSparams,trainExamplesOrgIndices,exampleInfo);
  
    %% select them
    [ntrainExamples,ntestExamples,featuresToKeep] = aux_selectFeatures(trainExamples,testExamples,nVoxelsToUse,sortedFeatures);
  else
    % no feature selection
    ntrainExamples = trainExamples;
    ntestExamples  = testExamples;
    sortedFeatures = 1:nFeatures;
  end
  nVoxelsToUse = size(ntrainExamples,2);
  nTestPerFold(k)   = size(testExamples,1); 
  
  %
  % outlier detection computation, computes score for each example
  % in the test set
  %
  
  switch method
    
   case {'SVD'}
    
    if subtractFeatureMean
      ntrainExamples = ntrainExamples - repmat(mean(ntrainExamples,1),size(ntrainExamples,1),1);
    end

    [u,s,v] = compute_fastSVD(ntrainExamples);

    d = diag(s);
    d = d .^ 2; dtotal = sum(d);
    
    tmp = cumsum(d/dtotal);
    tmp = find(tmp >= percentageOfVariance);
    n   = tmp(1); % # of components to keep

    u   = u(:,  1:n);
    s   = s(1:n,1:n);
    v   = v(:,  1:n);
    
    ntrainExamples = u*s; clear u;
    w = v';
    
    % now express the tested examples in terms of the v basis
    usTest = (inv(w*w' + speye(n)*1e-10) * w * ntestExamples')';clear w;
    
    rtestExamples = usTest*v';

    scoresPerFold{k} = sum( (rtestExamples-ntestExamples) .^ 2, 2);
   
   case {'nbayesPooled','logisticRegression2'}
  
    trainedClassifier   = trainClassifier(ntrainExamples,trainLabels,method,params);
    % store scores for fold/# of features combination
    tmpScores = applyClassifier(ntestExamples,trainedClassifier);
    
    scoresPerFold{k} = max(tmpScores,[],2);
%    scoresPerFold{k} = max(tmpScores,[],2) - min(tmpScores,[],2);

   case {'random'}
    scoresPerFold{k} = zeros(nTestPerFold(k),1);
    
   otherwise
    fprintf('%s: method %s is not supported\n',this,method);pause;return;
  end    
    
end; % end of loop over folds

%% for each number of features, assemble a score array

scores = zeros(nExamples,1);
eidx = 1;      
for k = 1:nFolds
  scores(eidx:(eidx+nTestPerFold(k)-1),:) = scoresPerFold{k};
  eidx = eidx+nTestPerFold(k);
end
    
%imagesc([examples,repmat(scores,1,10)]);

%
% Compute example outlier rankings (first is most outlying)
%

    
switch method
  
 case {'SVD'}
  
  [exampleScores,exampleRanking] = sort(scores);
  exampleRanking = flipud(exampleRanking); % large is outlier
  exampleScores  = flipud(exampleScores);
  
 case {'nbayesPooled','logisticRegression2'}
    
  [exampleScores,exampleRanking] = sort(scores);   % small is outlier
  
 case {'random'}
  % rank from random sort of the examples
  exampleRanking = randperm(nExamples);
  exampleScores  = exampleRanking;
  
 otherwise
  fprintf('%s: method %s is not supported\n',this,method);pause;return;
end




%% figure out which examples belong to which condition, and group
%% prior to dividing them into folds

% Divides the examples of each class into blocks with <exampleBlockSize> examples,
% which will end up together in cross-validation folds. By default,
% <exampleBlockSize> == 1, which corresponds to doing a balanced
% leave-1-out, with 1 example of each class in a fold
%
% <exampleBlockSize> == 0 is a special setting used to do unbalanced leave-1-out,
% in the presence of different numbers of examples of each
% condition. It will find out the condition with the fewest
% examples, and create 1 block for each of those.

function  [exampleBlocks] = aux_structureExamples(varargin)

DEBUG  = 0;
this   = 'aux_structureExamples';
labels = varargin{1};
if nargin; exampleBlockSize = varargin{2}; else; exampleBlockSize = 1; end

tmp = find(labels > 0);
conditions = unique(labels(tmp)); nConds = length(conditions);

%% find out how many examples there are of each cond and what indices they're at
for c = 1:nConds
  cond  = conditions(c);  
  examplesCond{c}  = find(labels == cond);
  nExamplesCond(c) = length(examplesCond{c});
end  

% HACK to have fewer examples in one of the classes
%labels(examplesCond{1}(1:3)) = 0;
%examplesCond{1}  = find(labels == conditions(1));
%nExamplesCond(1) = length(examplesCond{1});

minExamplesCond = min(nExamplesCond);
nFolds          = minExamplesCond;

%% now assign them to blocks
for c = 1:nConds
  
  if exampleBlockSize
    
    if rem(nExamplesCond(c),exampleBlockSize)
      fprintf('%s: error: #examples=%d in c=%d must be divisible by %d\n',this,nExamplesCond(c),cond,exampleBlockSize); pause; return
    else
      nBlocks(c) = nExamplesCond(c) / exampleBlockSize;
    end

    % rows in original matrix of examples that will end up in each block
    bidx = 1;
    exampleBlocks{c} = cell(nBlocks(c),1);
    for b = 1:nBlocks(c)
      exampleBlocks{c}{b} = examplesCond{c}(bidx:(bidx+exampleBlockSize-1));
      bidx = bidx + exampleBlockSize;
    end
    
  else
    % unbalanced leave 1 out with varying #s of examples

    nBlocks(c) = minExamplesCond;

    exampleBlocks{c}     = cell(nBlocks(c),1);
    exampleBlockSizeCond = floor(nExamplesCond(c)/nFolds);
    exampleBlockSizeRem  = rem(nExamplesCond(c),nFolds);

    bidx = 1;

    fprintf('c=%d\t%d\t%d\t%d\n',c,exampleBlockSizeCond,exampleBlockSizeRem,nBlocks(c));
    for b = 1:nBlocks(c)
      exampleBlocks{c}{b} = examplesCond{c}(bidx:(bidx+exampleBlockSizeCond-1))';
      bidx = bidx + exampleBlockSizeCond;
    end
    if exampleBlockSizeRem
      % last block gets remaining images
      exampleBlocks{c}{b} = [exampleBlocks{c}{b},examplesCond{c}(bidx:end)'];
    end
  end
end

if DEBUG
  fprintf('#examples/cond\t'); fprintf('%d ',nExamplesCond);fprintf('\n');
  fprintf('trial#\t');fprintf('%2d ',1:length(labels));fprintf('\n');
  fprintf('label \t');fprintf('%2d ',labels);fprintf('\n');
  for c = 1:nConds
    cond = conditions(c);
    fprintf('cond=%d #blocks=%d\n',cond,nBlocks(c));
    for b = 1:nBlocks(c)
      fprintf('\t');fprintf('%d\t',exampleBlocks{c}{b});
      fprintf('\t');fprintf('%2d ',labels(exampleBlocks{c}{b}));
      fprintf('\n');
    end
  end
  pause
end

% in balanced leave-1-out this is not supposed to happen
if 0
  if exampleBlockSize & sum(diff(nBlocks))
    fprintf('%s: task conditions have different numbers of trials\n',this);
    disp(conditions); disp(nBlocks); pause; return;
  end
end

%
% Create the division of examples into folds,  by determining
% what indices(rows) in the example matrix match the blocks of
% examples for each condition
% 

function [foldExamplesTrain,foldExamplesTest] = aux_divideIntoFolds(labels,crossValidationMethod,exampleBlockSize)

this = 'aux_divideIntoFolds';
DEBUG = 0;

tmp = find(labels > 0);
conditions = unique(labels(tmp)); nConds = length(conditions);

% Determine what indices(rows) in the example matrix make the blocks of examples for
% each condition, e.g. exampleBlocks{b} contains as many cells as blocks. Each
% cell has the indices of class <b> examples belonging to that block.
% 
% <exampleBlockSize> determines how many examples to put in a
% block. If doing a balanced leave-1-out, we'd have one per block,
% and one block of of each class in the test fold.
%
% <exampleBlockSize> == 0 is a special setting used to do balanced leave-1-out,
% in the presence of different numbers of examples of each
% condition. It will find out the condition with the fewest
% examples, and create 1 block for each of those.
exampleBlocks  = aux_structureExamples(labels,exampleBlockSize);

nBlocksPerCond = zeros(nConds,1);
for c = 1:nConds;
  nBlocksPerCond(c) = length(exampleBlocks{c});
end

% now divide the examples into folds; for each fold, create a set
% of indices for its training and test sets

switch crossValidationMethod
  
 case {'leave1outExample'}
  % leaves out a single block of examples, goes by condition
  % if exampleBlockSize=1, this is vanilla leave-1-out  
  nFolds = sum(nBlocksPerCond);
  
  fold = 1;
  for c = 1:nConds  
    for b = 1:nBlocksPerCond(c)
      
      % test examples
      foldExamplesTest{fold}  = exampleBlocks{c}{b};

      foldExamplesTrain{fold} = [];
      for c1 = 1:nConds  
	for b1 = 1:nBlocksPerCond(c1)
	  
	  if (c1 ~= c) | (b1 ~= b)
	    foldExamplesTrain{fold} = [foldExamplesTrain{fold},exampleBlocks{c1}{b1}];
	  end
	end
      end
      foldExamplesTrain{fold} = sort(foldExamplesTrain{fold});
      
      fold = fold + 1;
    end
  end
  
  
 case {'leave1out','leave1outUnbalanced'}

  % test or redo the exampleBlocks structure, as appropriate
  switch crossValidationMethod
   case {'leave1out'}
    if sum(diff(nBlocksPerCond))
      fprintf('%s: for %s we need a similar # of examples per class\n',this,crossValidationMethod);
      disp(conditions); disp(nBlocksPerCond); pause; return;
    end
    
   case {'leave1outUnbalanced'}
    % no conditions, this is actually being enforced by the example
    % structuring code. The only thing that can change over conditions
    % is the block size (and the last block can be larger, to catch
    % remains)
  end
  
  % balanced leave 1 out (leave out 1 example block (b) per class)
  
  for b = 1:nBlocksPerCond(c)
    
    foldExamplesTest{b} = [];      
    for c = 1:nConds
      foldExamplesTest{b} = [foldExamplesTest{b},exampleBlocks{c}{b}];
    end
    
    foldExamplesTrain{b} = [];
    for bt = 1:(b-1)
      for c = 1:nConds
	foldExamplesTrain{b} = [foldExamplesTrain{b},exampleBlocks{c}{bt}];
      end
    end
    for bt = (b+1):nBlocksPerCond
      for c = 1:nConds
	foldExamplesTrain{b} = [foldExamplesTrain{b},exampleBlocks{c}{bt}];
      end
    end    
    
  end
  nFolds = nBlocksPerCond;
        
  
 case {'evenOdd'}
  % same # of blocks per condition
  nBlocksPerCond = length(exampleBlocks{1});
  
  % interleaved 2 fold cross validation (odd examples of each
  % class are in one fold, even examples are in the other)       
  foldOdd  = []; foldEven = [];
  
  for b = 1:2:nBlocksPerCond
    for c = 1:nConds; foldOdd  = [foldOdd,exampleBlocks{c}{b}]; end
  end
  for b = 2:2:nBlocksPerCond
    for c = 1:nConds; foldEven = [foldEven,exampleBlocks{c}{b}]; end
  end
  
  foldExamplesTest{1}  = foldOdd;
  foldExamplesTrain{1} = foldEven;
  foldExamplesTest{2}  = foldEven;
  foldExamplesTrain{2} = foldOdd;    
  nFolds = 2;    

  
 case {'halfVShalf'}
  % same # of blocks per condition
  nBlocksPerCond = length(exampleBlocks{1});
  
  % interleaved 2 fold cross validation (odd examples of each
  % class are in one fold, even examples are in the other)       
  foldHalfA  = []; foldHalfB = [];

  halfBlock = floor(nBlocksPerCond/2);
  
  for b = 1:halfBlock
    for c = 1:nConds; foldHalfA = [foldHalfA,exampleBlocks{c}{b}]; end
  end
  for b = (halfBlock+1):nBlocksPerCond
    for c = 1:nConds; foldHalfB = [foldHalfB,exampleBlocks{c}{b}]; end
  end
  
  foldExamplesTest{1}  = foldHalfA;
  foldExamplesTrain{1} = foldHalfB;
  foldExamplesTest{2}  = foldHalfB;
  foldExamplesTrain{2} = foldHalfA;    
  nFolds = 2;    
  
end; % end of switch over folds

if DEBUG
  for k = 1:nFolds
    fprintf('fold %d\n',k);
    fprintf('test size = %d\ttrain size = %d\n',length(foldExamplesTest{k}),length(foldExamplesTrain{k}));
    fprintf('test  idx\t');   fprintf('%2d ',foldExamplesTest{k});fprintf('\n');
    fprintf('test  label\t'); fprintf('%2d ',labels(foldExamplesTest{k}));fprintf('\n');
    fprintf('train idx\t');   fprintf('%2d ',foldExamplesTrain{k});fprintf('\n');
    fprintf('train label\t'); fprintf('%2d ',labels(foldExamplesTrain{k}));fprintf('\n');
  end
  pause
end


%% feature ranking

function [sortedFeatures,sortedMeasure] = aux_rankFeatures( varargin )

l = length(varargin);
% example arguments
trainExamples = varargin{1};
trainLabels   = varargin{2};
[nExamples,nFeatures] = size(trainExamples);
conditions    = unique(trainLabels); nConds = length(conditions);

% feature selection arguments
method = varargin{3};
methodParameters   = varargin{4};
examplesOrgIndices = varargin{5};
exampleInfo        = varargin{6};

% other arguments
%equalizeCycles = varargin{5};
%featureToGroup = varargin{9};
equalizeCycles  = 0;

%% Rank features

%% rank features by one of the methods
switch method

 case {'CCBI'}
  % select features using one of the CCBI methods
  [sortedFeatures,sortedMeasure] = summarizeELE_rankByCCBImethod(trainExamples,trainLabels,methodParameters{1});
%  [sortedFeatures,sortedMeasure] = summarizeEXA_rankByCCBImethod(trainExamples,trainLabels,methodParameters{1},examplesOrgIndices,exampleInfo);  

  
 case {'nestedCV'}
  classifierToNest   = methodParameters{1};
  errorMeasureToNest = methodParameters{2};

  % GNB discriminability of each feature in CV within the training set
%  [sortedFeatures,sortedMeasure] = summarizeELE_rankByNestedCV(trainExamples,trainLabels,classifierToNest,errorMeasureToNest);
  [sortedFeatures,sortedMeasure] = summarizeEXA_rankByNestedCV(trainExamples,trainLabels,classifierToNest,errorMeasureToNest);
  
 case {'nestedCVpairwise'}
  
  % same thing as nestedCV, but considers discriminability for each
  % pair of classes, ranks features for each pair, and then picks
  % the 1st position of each ranking, the 2nd, etc
  %  [sortedFeatures,sortedMeasure] = summarizeELE_rankByPairwise(trainExamples,trainLabels);
[sortedFeatures,sortedMeasure] = summarizeEXA_rankByPairwise(trainExamples,trainLabels);
  
 case {'activityTtest'}

  % the test type to perform should be "mean > 0" for normal voxels
  % but "mean != 0" for voxels after "zero average cycle image" norm
  if equalizeCycles; testType = 'differentFromZero'; else; testType = 'greaterThanZero';end

  %  [sortedFeatures,sortedMeasure] = summarizeELE_rankByMeanTest(trainExamples,trainLabels,testType);
  [sortedFeatures,sortedMeasure] = summarizeEXA_rankByMeanTest(trainExamples,trainLabels,testType);
  
 case {'contribution'}

  % doing CV on the training set, see how much each feature
  % contributes to making each example correlate with the ones of
  % each class in the training set, optionally trading that off
  % with the average contribution to correlation with examples of
  % the wrong class
  
  if 1
    measure             = methodParameters{1};
    combination         = methodParameters{2};
    avgTrainingExamples = methodParameters{3};
    
    [sortedFeatures,sortedMeasure] = summarizeELE_rankByDistContrib(trainExamples,trainLabels,measure,combination,avgTrainingExamples);
  else
    % EXPERIMENTAL, DOESN'T WORK YET
    % build a command line from parameters - ultimately, all methods
    % should be coded more or less like this
    cmd = '[sortedFeatures,sortedMeasure] = summarizeELE_rankByDistContrib(trainExamples,trainLabels';
    if ~isempty(methodParameters)
      for mp = 1:length(methodParameters)
	cmd = [cmd,',''',methodParameters{mp},''''];      
      end
    end
    cmd = [cmd,');'];
    eval(cmd);
  end
    
 case {'contribToWCsimilarity'}
  
  similarityMeasure = methodParameters{1};
  combinationMethod = methodParameters{2};
  
  [sortedFeatures,sortedMeasure] = summarizeELE_rankByWithinClassContrib(trainExamples,trainLabels,similarityMeasure,combinationMethod);
 
 case {'random'}
  [sortedFeatures,sortedMeasure] = summarizeELE_rankByRandom(trainExamples,trainLabels);

end; % switch over methods


%% average examples of the same class

function [newExamples,newLabels] = aux_averageSameLabelExamples(examples,labels)

conditions = unique(labels); nConds = length(conditions);
nLabels = size(labels,2); nFeatures = size(examples,2);

newExamples = zeros(nConds,nFeatures);
newLabels   = zeros(nConds,nLabels);

for c = 1:nConds
  cond  = conditions(c);  
  examplesCond{c}  = find(labels == cond);
  nExamplesCond(c) = length(examplesCond{c});

  newExamples(c,:) = mean(examples(examplesCond{c},:),1);
  newLabels(c,1)   = conditions(c);
end  


%% select features on both train and test set according to the ranking

function [trainExamples,testExamples,featuresToKeep] = aux_selectFeatures(trainExamples,testExamples,nToKeep,sortedFeatures)

nFeatures      = size(trainExamples,2);
if nToKeep > nFeatures; nToKeep = nFeatures; end
featuresToKeep = sortedFeatures(1:nToKeep);
trainExamples  = trainExamples(:,featuresToKeep);
testExamples   = testExamples(:,featuresToKeep);


%% figure out which examples belong to which condition, and group
%% prior to dividing them into folds

% Divides the examples of each class into blocks with <exampleBlockSize> examples,
% which will end up together in cross-validation folds. By default,
% <exampleBlockSize> == 1, which corresponds to doing a balanced
% leave-1-out, with 1 example of each class in a fold
%
% <exampleBlockSize> == 0 is a special setting used to do unbalanced leave-1-out,
% in the presence of different numbers of examples of each
% condition. It will find out the condition with the fewest
% examples, and create 1 block for each of those.


function [] = test()

% Create test examples for this. All the same except
% each epoch has a single outlier. In addition, half
% the features are noise.
% 

nFeatures = 100;
nClasses  = 14;
nEpochs   = 6;

labels   = repmat((1:nClasses)',nEpochs,1);
examples = repmat(randn(1,nFeatures),nClasses*nEpochs,1);

outliers = [];
idx      = 1;
for e = 1:nEpochs
  examples(idx,:) = randn(1,nFeatures);
  outliers = [outliers,idx];
  idx = idx + nClasses + 1;
end
nExamples = size(examples,1);

outliers
examples = [examples,randn(nExamples,nFeatures/2)];

