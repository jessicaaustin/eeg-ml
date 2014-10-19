% applyClassifier(examples,classifier)
%
% Apply a learnt classifier (from trainClassifier) to the new examples and
% get the scores for classification.
%
% This function outputs a matrix of scores of the labelling over all the
% classes present in the test set. The columns of that matrix
% correspond to the labels sorted in increasing order of their
% numeric values
% (e.g. columns 1-6 would be labels 2-7 in the sixcategories study,
% which has 6 classes)
%
% The specific scoring measure will depend on the classifier used,
% see classifier code for details
%
% Input:
%  - examples, a matrix of real numbers with dimenstion of #examples*#features
%  - classifier (learnt by trainClassifier)
%
% Output:
%  - score, a matrix of real number with dimension of #examples*#classes,
%    Each row contains [score(1st class|data) ... score(kth class|data)]
%
% Dependencies:
%  - serious dependence on the "classifier" format produced by trainClassifier.
%
% Example:
%  - [predictions] = applyClassifier(examples,classifier);
%
% History:
%  - Oct07,2005 Wei - redocument
%  - 10 Mar 05 - fpereira - adapted from previous code



function [scores] = applyClassifier( varargin )

%
% Process parameters
% 

l = length(varargin);
if l < 2; help applyClassifier; return; end;
  
examples          = varargin{1};
trainedClassifier = varargin{2};

classifierType   = trainedClassifier.classifier;
classifierParams = trainedClassifier.classifierParameters;
models           = trainedClassifier.models;

lastpos           = size(models,1);
trainingSetInfo   = trainedClassifier.trainingSetInfo;
sortedLabelValues = trainingSetInfo.sortedLabelValues;
nClasses          = length(sortedLabelValues);
nFeatures         = size(examples,2);
nExamples         = size(examples,1);
  
% This matrix will keep the prediction scores and feature relevance values
% corresponding to each example
% - scores is a #examples*#classes matrix
%   (each line is scores for one example, in label value order)
% - cFeatureContribution is a cell array where each cell contains
%   a featureContribution map for a class
%   A featureContribution map is a #examples*#features matrix
  
scores            = zeros(nExamples,nClasses);
cFeatureContribution = cell(nClasses,1);
cFeatureLogProbs  = cell(nClasses,1);
    
%
% run the classifier - eventually we will encapsulate this and
% the code for different kinds of classifier, for now a large switch statement

trainingSetInfo = models{nClasses+2};
priors          = trainingSetInfo.classPriors;
  
switch classifierType

  %
  % Logistic regression
  %
  
  
 case {'logisticRegression2'}

  % discriminative model weights;
  W = models{nClasses+1};
  
  tmp = models{nClasses+3}{1}; % parameters used to train
  transform = tmp{3}; lambda = tmp{2}; eta = tmp{1};

  switch transform
    
   case {'none'}
    % vanilla logistic regression
    tmp    = exp([ones(nExamples,1),examples] * W);
    scores = tmp ./ repmat(sum(tmp,2),1,size(W,2));

    
   case {'correlation','hybrid'}
    
    %% normalize each example by subtracting its mean and / by its standard deviation
    eachExampleMean = mean(examples,2);
    eachExampleStdv = std(examples,0,2);  
    examples = examples -  repmat(eachExampleMean,1,nFeatures);
    examples = examples ./ repmat(eachExampleStdv,1,nFeatures);
      
    %% get mean examples from the training set (one per class)
    %% (already normalized)
    meanExamples = models{nClasses+3}{2};
    
    %% produce per class score
    for c = 1:nClasses
      tmp = [ones(nExamples,1),(examples .* repmat(meanExamples(c,:),nExamples,1))];
      scores(:,c) = tmp * W(:,c);
      %      scores(:,c) = sum(tmp,2);
    end
    scores = scores ./ repmat(sum(scores,2),1,nClasses); % normalize

    
  end
    
  
  %
  % Logistic regression
  %

 case {'logisticRegression', 'SMLR'}
  weights=models{lastpos-2};
  examples=[ones(size(examples,1),1) examples];
  
  tmp = exp(examples * weights);
  scores = tmp ./ repmat(sum(tmp,2),1,size(weights,2));
  
  % <Lucas Tan>
  % As the product of examples*weights might produce a large value, 
  % say 1000, an element in tmp could be e^1000 = Infinity (overflow).
  % Recall the logistic function is e^(wx) / (1 + e^(wx)),
  % And so the corresponding scores(i) would be Inf/(1+Inf) == NaN.
  % To fix it, we replace all NaN's with the correct value of 1.
  % The probability of the other class would be 1/(1+Inf) ==0 
  % which would be correctly handled by Matlab.

  global DEBUG;
  if exist('DEBUG', 'var') && ~isempty(strfind(DEBUG, '-LRFINPREC')),
	scores(isnan(scores)) = 1;
  end
  % </Lucas Tan>

  switch classifierType
  case { 'L1logisticRegression','L2logisticRegression'}
      % swap the columns to ensure correct labeling
      scores=[scores(:,2) scores(:,1)];
  end
  
  %
  % Gaussian naive bayes and friends
  %
  
 case {'nbayes','nbayesPooled','nbayes-unitvariance'}
    
  % Compute log(p(class c|example)) for each example
  % and create the score table for this method
  % (a #examples*#classes matrix where higher score means higher
  % likelihood that that is the correct class)

  featureWeights  = ones(1,nFeatures); % for now
  
  for c=1:nClasses
    
    classMeans = models{c}{1};
    classStdvs = models{c}{2};
    classPrior = priors(c);
    
    % compute ln(p(class|data))
    meanMask = repmat(classMeans,[nExamples,1]);
    stdvMask = repmat(classStdvs,[nExamples,1]);
    
    % compute per feature p(feature|class), for all features and
    % examples at once
    tmp = examples;
    tmp = tmp -  meanMask;
    tmp = tmp ./ stdvMask;
    tmp = -1 * (tmp .^ 2);
    tmp = tmp / 2;
    tmp = tmp - log(sqrt(2*pi)*stdvMask);
    
    % weigh each feature contribution and add them up
    scores(:,c) = tmp * featureWeights' + log(classPrior);
  end
  
  %
  % k nearest neighbours
  %
  
 case {'knn'}
    
  trainingSetInfo = models{nClasses+2};
  classPriors     = trainingSetInfo.classPriors;
  logClassConst   = zeros(size(classPriors));  
  
  % the model in this case is the training data
  % models{nClasses+1}{1}; % examples
  % models{nClasses+1}{2}; % labels
  % models{nClasses+1}{3}; % k 
  k = models{nClasses+1}{3}; % k 
  distance = models{nClasses+1}{4}; % distance
  
  scores = classifierKNN2(k,models{nClasses+1}{1},models{nClasses+1}{2},examples,distance);
    
  %
  % SVM
  %
  
 case {'svmlight'}
  % uses the MATLAB wrapper for Thorsten Joachim's SVMlight
  % See SVM/svml.m for details on options available. For now,
  % we'll have to settle for picking a kernel and its parameters
  % (follows SVMdemsvml1.m example)
  
  trainingSetInfo   = models{nClasses+2};
  classPriors       = trainingSetInfo.classPriors;
  nClasses           = trainingSetInfo.nClasses;
  sortedLabelValues = trainingSetInfo.sortedLabelValues;
  
  if nClasses > 2
    fprintf('ERROR: current SVM code only allows 2 class problems or a pairwise voting classifier based on SVM\n');
    pause; return;
  end
    
  ypred = svmlfwd( models{nClasses+1}, examples );
  
  % each entry in ypred is negative (-1 prediction) or positive
  % (+1 prediction) , convert to 1 2
  indicesA = find(ypred <= 0);
  indicesB = find(ypred > 0);
  ypred(indicesA) = 1;
  ypred(indicesB) = 2;
    
  % Create a score array with one row per example, with one
  % column per label (in sortedLabelValues order). In that row,
  % the predicted label has value 1 and the others 0
  scores = zeros(nExamples,nClasses);
  indices = sub2ind([nExamples,nClasses],(1:nExamples)',ypred);
  scores(indices) = 1;
  
    
  %
  % Neural networks
  %
  
 case {'neural'}
  % uses the NETLAB code by Ian Nabney 
  % See Netlab/net.m for details on options available
  
  % get trained network and information
  trainingSetInfo   = models{nClasses+2};
  classPriors       = trainingSetInfo.classPriors;
  nClasses          = trainingSetInfo.nClasses;
  sortedLabelValues = trainingSetInfo.sortedLabelValues;
  
  net      = models{nClasses+1}; % trained network
  examples = normalize(examples); % normalize the data to mean 0 std 1
  netType  = net.ourNetType; % HACK to pass the network type with it
  
  % apply - scores are posterior prob of each class given 
  fprintf('applying a %s neural network\n',netType);
  switch netType
   case {'MLP'}
    scores = mlpfwd(net,examples);
   case {'GLM'}
    scores = glmfwd(net,examples); 
  end
  

  %
  % LDA/QDA, plus variations involving preprocessing with SVD or CV
  %
  
 case {'lda','qda','ldaSVD','ldaCV','qdaSVD','qdaCV'}
  
  covarianceMatrix = models{nClasses+3}{1};
    
  switch classifierType
    
   case {'ldaSVD','qdaSVD'}
    % V matrix from SVD projects data into lowdim subspace
    projectionMatrix = models{nClasses+3}{2};
    examples  = examples * projectionMatrix;
    nFeatures = size(examples,2); 
    
   case {'ldaCV','qdaCV'}
    % V matrix from initial SVD
    projectionMatrixSVD = models{nClasses+3}{2}{1};
    examples  = examples * projectionMatrixSVD;
    % V matrix from subsequent CV
    projectionMatrixCV  = models{nClasses+3}{2}{2};
    examples  = examples * projectionMatrixCV;
    nFeatures = size(examples,2); 
    
  end
    
  for c = 1:nClasses
    
    classMeans = models{c}{1};
    classPrior = priors(c);
    
    % log( p(example | class) )
    switch classifierType
     case {'lda','ldaSVD','ldaCV'}
      %	fprintf('%s %d %d\n',classifierType,size(examples));
      size(covarianceMatrix)
      rank(covarianceMatrix)
      rcond(covarianceMatrix)
      tmp = mymvnpdf( examples, classMeans, covarianceMatrix );
      
     case {'qda','qdaSVD','qdaCV'}
      % each class has a different covariance matrix
      %	fprintf('%s %d %d\n',classifierType,size(examples));
      tmp = mymvnpdf( examples, classMeans, covarianceMatrix{c} );
    end
    
    tmp = log( tmp );
    tmp = tmp + log(classPrior);
    
    scores(:,c) = tmp;
  end

  
  %
  % pairwise voting (voting scheme described below)
  %
        
 case {'pairwise'}
  
  trainingSetInfo           = models{nClasses+2};
  classifierParameters      = trainingSetInfo.classifierParameters;
  nParameters               = length(classifierParameters);
  classifierToUse           = classifierParameters{1}
  classifierParametersToUse = {};
  if length(classifierParametersToUse) > 1
    classifierParametersToUse = classifierParameters{2};
  end
  
  allPairModels = models{nClasses+1};
  allPairScores = zeros(nExamples,nClasses);
  
  fprintf('pairwise classifier with %s\n',classifierToUse);
  
  % 1) Apply each pairwise classifier on all the examples, yielding labels P1/P2
  % 2) For each example, add 1 to the vote count for either P1 or P2
  % 3) At the end, label according to the label with more counts    
  
  for l1 = 1:nClasses-1
    for l2 = l1+1:nClasses
      
      c1 = sortedLabelValues(l1);
      c2 = sortedLabelValues(l2);
      labelPair = [c1 c2];
      indexPair = [l1 l2];
      fprintf('\t%d %d\n',c1,c2);
      
      % apply it -> scores over that condition pair
      scores = applyClassifier(examples,allPairModels{l1,l2},classifierToUse,classifierParametersToUse);
      % make a classification decision based on picking, for each
      % example, the label in the pair that has a higher score
      [sortedScores,sortedIndices] = sort(scores,2);
      highScore = sortedIndices(:,2);
      
      % update the vote count of the picked class for each example
      highIdx   = indexPair(highScore);
      
      % not efficient, but this has to be done quickly
      for e = 1:nExamples
	allPairScores(e,highIdx(e)) = allPairScores(e,highIdx(e)) + 1;
      end
      
    end; % label 1
  end; % label 2
  
  % normalize to % of choices
  voteTotals = sum(allPairScores,2);
  scores     = allPairScores ./ repmat(voteTotals,1,nClasses);
  

  
 case {'nnets'}
  scores = sim(models{nClasses+1},examples');
  % now normalize the scores
  size(scores);
  for i=1:1:size(examples,1)
    scores(:,i) = scores(:,i) / sum(scores(:,i));
  end

  
  
 case {'svm'}
  
  trainingSetInfo = models{nClasses+2};
  classPriors     = trainingSetInfo.classPriors;
  logClassConst   = zeros(size(classPriors));  
  
  % the model in this case is the training data
  % models{nClasses+1}{1}; % examples
  % models{nClasses+1}{2}; % labels
  % models{numLabels+1}{3} = ker; % kernel type
  % models{numLabels+1}{4} = alpha;
  % models{numLabels+1}{5} = bias;
  % arguments are
  scores = classifierSVM2(models{nClasses+1}{1},models{nClasses+1}{2},examples,models{nClasses+1}{3},models{nClasses+1}{4},models{nClasses+1}{5});
  

  %
  % EXPERIMENTAL CODE, LEAVE ALONE AND DON'T LOOK TOO CLOSELY
  % (fpereira)
  
  
 case {'nbayesPooledResampling'}
  
  % compute resampled feature weights, which are *by example*
  % training set is models{nClasses+3}{1};
  % training labels is models{nClasses+3}{2};
  [exampleFeatureWeights] = compute_testSetRatio(models{nClasses+3}{1},models{nClasses+3}{2},examples);
  
  % Compute log(p(class c|example)) for each example
  % and create the score table for this method
  % (a #examples*#classes matrix where higher score means higher
  % likelihood that that is the correct class)
  
  for c=1:nClasses
    
    classMeans = models{c}{1};
    classStdvs = models{c}{2};
    classPrior = priors(c);
    
    % compute ln(p(class|data))
    meanMask = repmat(classMeans,[nExamples,1]);
    stdvMask = repmat(classStdvs,[nExamples,1]);
    
    tmp = examples;
    tmp = tmp -  meanMask;
    tmp = tmp ./ stdvMask;
    tmp = -1 * (tmp .^ 2);
    tmp = tmp / 2;
    tmp = tmp - log(sqrt(2*pi)*stdvMask);
    
    % weigh each feature contribution and add them up
    scores(:,c) = sum(tmp .* exampleFeatureWeights,2) + log(classPrior);
  end
    

end
  
%
% Ancillary functions or functions that Matlab R13 doesn't have yet
%
  
  
%% Normalize each feature to have mean 0 and standard deviation 1

function [Y] = normalize(X)

[nExamples,nFeatures] = size(X);
meanX = mean(X,1);
stdvX = std(X,0,1);

Y = X -  repmat(meanX,[nExamples,1]);
Y = Y ./ repmat(stdvX,[nExamples,1]);

%
% Copy of the MATLAB function so that versions earlier than R14 can
% run this. REPLACE all mentions of this function by the official
% mvnpdf as soon as we're up to R14.
% 

function y = mymvnpdf(X, Mu, Sigma)

if nargin < 1 | isempty(X)
    error('stats:mymvnpdf:InputSizeMismatch','Requires the input argument X.');
elseif ndims(X) > 2
    error('stats:mymvnpdf:InvalidData','X must be a matrix.');
end

% Get size of data.  Column vectors provisionally interpreted as multiple scalar data.
[n,d] = size(X);

% Assume zero mean, data are already centered
if nargin < 2 | isempty(Mu)
    X0 = X;

% Get scalar mean, and use it to center data
elseif prod(size(Mu)) == 1
    X0 = X - Mu;

% Get vector mean, and use it to center data
elseif ndims(Mu) == 2
    [n2,d2] = size(Mu);
    if d2 ~= d % has to have same number of coords as X
        error('stats:mymvnpdf:InputSizeMismatch',...
              'X and MU must have the same number of columns.');
    elseif n2 == n % lengths match
        X0 = X - Mu;
    elseif n2 == 1 % mean is a single row, rep it out to match data
        X0 = X - repmat(Mu,n,1);
    elseif n == 1 % data is a single row, rep it out to match mean
        n = n2;
        X0 = repmat(X,n2,1) - Mu;
    else % sizes don't match
        error('stats:mymvnpdf:InputSizeMismatch',...
              'X or MU must be a row vector, or X and MU must have the same number of rows.');
    end
    
else
    error('stats:mymvnpdf:BadMu','MU must be a matrix.');
end

% Assume identity covariance, data are already standardized
if nargin < 3 | isempty(Sigma)
    % Special case: if Sigma isn't supplied, then interpret X
    % and Mu as row vectors if they were both column vectors
    if d == 1 & prod(size(X)) > 1
        X0 = X0';
        [n,d] = size(X0);
    end
    xRinv = X0;
    sqrtInvDetSigma = 1;
    
% Single covariance matrix
elseif ndims(Sigma) == 2
    % Special case: if Sigma is supplied, then use it to try to interpret
    % X and Mu as row vectors if they were both column vectors.
    if (d == 1 & prod(size(X)) > 1) & size(Sigma,1) == n
        X0 = X0';
        [n,d] = size(X0);
    end
    
    % Make sure Sigma is the right size
    if size(Sigma,1) ~= d | size(Sigma,2) ~= d
        error('stats:mymvnpdf:BadSigma',...
              'SIGMA must be a square matrix with size equal to the number of columns in X.');
    else
        % Make sure Sigma is a valid covariance matrix
        [spd,R] = myisspd(Sigma);
        if spd
            % Create array of standardized data, vector of inverse det
            xRinv = X0 / R;
            sqrtInvDetSigma = 1 / prod(diag(R));
        else
            error('stats:mymvnpdf:BadSigma',...
                  'SIGMA must be symmetric and positive definite.');
        end
    end
    
% Multiple covariance matrices
elseif ndims(Sigma) == 3
    % Special case: if Sigma is supplied, then use it to try to interpret
    % X and Mu as row vectors if they were both column vectors.
    if (d == 1 & prod(size(X)) > 1) & size(Sigma,1) == n
        X0 = X0';
        [n,d] = size(X0);
    end
    
    % Data and mean are a single row, rep them out to match covariance
    if n == 1 % already know size(Sigma,3) > 1
        n = size(Sigma,3);
        X0 = repmat(X0,n,1); % rep centered data out to match cov
    end

    % Make sure Sigma is the right size
    if size(Sigma,1) ~= d | size(Sigma,2) ~= d
        error('stats:mymvnpdf:InputSizeMismatch',...
              'Each page of SIGMA must be a square matrix with size equal to the number of columns in X.');
    elseif size(Sigma,3) ~= n
        error('stats:mymvnpdf:InputSizeMismatch',...
              'SIGMA must have one page for each row of X.');
    else
        
        % Create array of standardized data, vector of inverse det
        xRinv = zeros(n,d);
        sqrtInvDetSigma = zeros(n,1);
        for i = 1:n
            % Make sure Sigma is a valid covariance matrix
            [spd,R] = myisspd(Sigma(:,:,i));
            if spd
                xRinv(i,:) = X0(i,:) / R;
                sqrtInvDetSigma(i) = 1 / prod(diag(R));
            else
                error('stats:mymvnpdf:BadSigma',...
                      'SIGMA must be symmetric and positive definite.');
            end
        end
    end
   
elseif ndims(Sigma) > 3
    error('stats:mymvnpdf:BadSigma',...
          'SIGMA must be a matrix or a 3 dimensional array.');
end

% Exponents in pdf are the inner products of the standardized data
quadform = sum(xRinv.^2, 2);
y = sqrt((2*pi)^(-d)) * sqrtInvDetSigma .* exp(-0.5*quadform);


%
% More matlab R13/R14 functions
%

function [t,R] = myisspd(Sigma)

% Test for square, symmetric
[n,m] = size(Sigma);
if (n == m) & all(all(abs(Sigma - Sigma') < 10*myeps(max(abs(diag(Sigma))))))
    % Test for positive definiteness
    [R,p] = chol(Sigma);
    if p == 0
        t = 1;
    else
        t = 0;
    end
else
    R = [];
    t = 0;
end


function [r] = myeps(x)

E = floor(log2(x));
if isa(x,'single')
     r = 2^(E-23);
elseif isa(x,'double')
     r = 2^(E-52);
else
error('myeps: argument must be a single or a double');
end
