% trainClassifier
%
% Trains a classifier on a set of examples
%
% This code wraps around the code for several different classifiers
% and knows how to package their results in a uniform manner, to
% pass to applyClassifier.m
%
% Input:
% - examples (#examples x #features matrix)
% - labels   (#examples x 1 vector)
% - classifier (classifier name, e.g. 'nbayes')
% - parameters, varied for different classifiers
%
% Output:
% - a "trainedClassifier" structure
%
%
% models = cell(nClasses+3,1);
%
% Cells 1 to nclasses - contain the generative model for each class
% Cell nclasses+1 - contains the discriminative model, if any
% Cell nclasses+2 - contains information about training set
% Cell nclasses+3 - contains any extra information, model specific
%
% Generative model - each cell contains a cell array
% where each cell has one parameter, e.g.
%
%  models{1}{1} = mean class 1
%  models{1}{2} = stdev class 1
%  models{1}{3} = prior prob class 1
%
% Discriminative model - a cell array of sets of weights
%  models{nclasses+1} = cell(1,1);
%  
% Information
%  models{nclasses+2} = cell(1,1);
% containing
%  trainingSetInfo.nExamples (for the time being)
%
% Other classifier specific information
%  models{nclasses+3}
%
% Dependencies:
%
% Example:
%   classifier = trainClassifier(examples, labels, 'logisticRegression', {0.01 0.001 10});
%
% History
% - Oct 07,2005 Wei - redocument
% - 10 Mar 05 - fpereira - adapted from several past versions
% - 01 Jun 05 - wew - add logisticRegression

function [trainedClassifier] = trainClassifier( varargin )

%
% process examples
%

l = length(varargin);
if l < 3; help trainClassifier; return; end
  
examples       = varargin{1};
labels         = varargin{2};
classifierType = varargin{3};

classifierParameters = {};
if l > 3; classifierParameters = varargin{4}; end
  
% figure out a few things
sortedLabelValues      = sort(unique(labels));
nClasses              = length(sortedLabelValues);
[numTrain,numFeatures] = size(examples);

% find the indices of examples with each label
for l = 1:nClasses
  label = sortedLabelValues(l);
  examplesWithLabel{l} = find(labels == label);
end


%
% Train the classifier
%

models = cell(nClasses+3,1);
lcp = length(classifierParameters);

switch classifierType
  
 case {'logisticRegression'}
  %fprintf('trainClassifier: using %s with parameters\n',classifierType);disp(classifierParameters);  
  [models] = classifierLogisticRegression(examples,labels,classifierParameters);
  
 case {'SMLR'}
  %fprintf('trainClassifier: using %s with parameters\n',classifierType);disp(classifierParameters);  
  [models] = classifierSMLR(examples,labels,classifierParameters);
  
 case {'nbayes','nbayesPooled','lda','qda','ldaSVD','ldaCV','qdaSVD','qdaCV','nbayes-unitvariance','nbayesPooledResampling'}
    
  %fprintf('trainClassifier: using %s with parameters\n',classifierType);disp(classifierParameters);  
  [models] = classifierBayes(examples,labels,classifierType,classifierParameters);

  
 case {'knn'}
    
  % defaults
  k = 1;
  distance = 'euclidean';

  if lcp > 0
    % override defaults
    k = classifierParameters{1};
    if lcp > 1
      distance = classifierParameters{2};
    end
  end
  
  %fprintf('trainClassifier: knn: using k=%d and %s distance\n',k,distance);

  % store the models in the "discriminative" part
  models{nClasses+1} = cell(3,1);
  models{nClasses+1}{1} = examples;
  models{nClasses+1}{2} = labels;
  models{nClasses+1}{3} = k;
  models{nClasses+1}{4} = distance;

    
 case {'svmlight'}
  
  % uses the MATLAB wrapper for Thorsten Joachim's SVMlight
  % See SVM/svml.m for details on options available. For now,
  % we'll have to settle for picking a kernel and its parameters
  % (follows SVM/demsvml1.m example)
  
  %   'Kernel'         -t       {0..4}, default value 1
  %                             Type of kernel function:
  %                             0: linear
  %                             1: polynomial (s a*b+c)^d
  %                             2: radial basis function exp(-gamma ||a-b||^2)
  %                             3: sigmoid tanh(s a*b + c)
  %                             4: user defined kernel from kernel.h
  %   'KernelParam'    -d, -g, -s, -r, -u
  %                             Depending on the kernel, this vector
  %                             contains [d] for polynomial kernel, [gamma]
  %                             for RBF, [s, c] for tanh kernel, string for
  %                             user-defined kernel
  
  %fprintf('trainClassifier: using %s with parameters\n',classifierType);disp(classifierParameters);
  
  if lcp
    kernel       = classifierParameters{1}; % kernel type
    kernelParams = classifierParameters{2}; % a vector of parameter values    
  else
    % default to a linear kernel
    kernel = 0;
    kernelParams = [];
  end
  net = svml('tmpsvm', 'Kernel', kernel , 'KernelParam', kernelParams);
    
  % convert labels to +1 or -1
  if nClasses > 2
    fprintf('ERROR: current SVM code only allows 2 class problems or a pairwise voting classifier based on SVM\n');
    pause
    return;
  else
    indices1 = find(labels == sortedLabelValues(1));
    indices2 = find(labels == sortedLabelValues(2));
    labels(indices1) = -1;
    labels(indices2) = 1;
  end
    
  % store model in the "discriminative" part
  models{nClasses+1} = svmltrain(net, examples, labels);
    
    
 case {'neural'}
  
  % uses the NETLAB code by Ian Nabney 
  % See Netlab/net.m for details on options available
  net = classifierNeuralNetwork(examples,labels,classifierParameters);

  % store model in the "discriminative" part
  models{nClasses+1} = net;


 case {'pairwise'}
  
  % Special classifier that trains a given classifier on
  % every pair of classes and then uses all those models
  % to reach a decision by voting
  classifierToUse           = classifierParameters{1};
  classifierParametersToUse = {};
  if length(classifierParametersToUse) > 1
    classifierParametersToUse = classifierParameters{2};
  end
  
  allPairModels             = cell(nClasses,nClasses);
  
  fprintf('pairwise classifier with %s\n',classifierToUse);
  
  for l1 = 1:nClasses-1
    for l2 = l1+1:nClasses
      
      c1 = sortedLabelValues(l1);
      c2 = sortedLabelValues(l2);
      fprintf('\ttraining %d %d\n',c1,c2);
      
      % separate the data of those two conditions
      indices      = find((labels==c1)|(labels==c2));
      pairExamples = examples(indices,:);
      pairLabels   = labels(indices,:);
      
      % train a classifier on data of the two conditions
      allPairModels{l1,l2} = trainClassifier(pairExamples,pairLabels,classifierToUse,classifierParametersToUse);
    end; % label 1
  end; % label 2
    
  % store the cell array with models in the "discriminative" part
  models{nClasses+1} = allPairModels;

  
 case {'nnets'}
  [models] = classifierNNets(examples,labels);

  
 case {'svm'}
  % the old SVM code (inefficient, implementation for a class)

  % other Parameters for SVMs
  params = size(classifierParameters,2);
    
  % default ker
  ker = 'linear';
  if (params > 0); ker = classifierParameters{1};
  end;
  % default C - bound on multipliers
  C = inf; 
  if (params > 1); C = classifierParameters{2}; end;
  %default p1
  global p1;
  p1 = 1;
  if (params > 2); p1 = classifierParameters{3};end;
  % default p2 
  global p2;
  p2 = 0;
  if (params > 3); p2 = classifierParameters{4}; end;
  
  % convert the labels to -1 and 1
  % WARNING: assumes there are only two classes
  correctedTrainLabels = zeros(size(labels));
  indices = find(labels == sortedLabelValues(1));
  correctedTrainLabels(indices) = -1;
  indices = find(labels == sortedLabelValues(2));
  correctedTrainLabels(indices) = 1;
  
  % train it
  [nsv,alpha,bias] = svc(examples,correctedTrainLabels,ker,C);
  
  % store the models in the "discriminative" part    
  models = cell(nClasses+3,1);
  models{nClasses+1} = cell(5,1);
  models{nClasses+1}{1} = examples;
  models{nClasses+1}{2} = labels;
  models{nClasses+1}{3} = ker;
  models{nClasses+1}{4} = alpha;
  models{nClasses+1}{5} = bias;
  
 otherwise
  fprintf('trainClassifier: error: classifier %s is not supported\n',classifierType);
end
  
%
% Store training set information
%

switch classifierType
 case {'knn','svmlight','pairwise','neural','nnets','svm','logisticRegression'}

  % Training Set information
  trainingSetInfo.classifierParameters = classifierParameters;
  trainingSetInfo.nExamples            = numTrain;
  trainingSetInfo.nFeatures            = numFeatures;
  trainingSetInfo.nClasses             = nClasses;
  trainingSetInfo.sortedLabelValues    = sortedLabelValues;
  trainingSetInfo.classPriors          = zeros(nClasses,1);
  for l=1:nClasses
    trainingSetInfo.classPriors(l) = length(find(labels==sortedLabelValues(l)));
  end
  trainingSetInfo.classPriors = trainingSetInfo.classPriors/numTrain;            
  
 otherwise
  % the bayesian classifiers set it up inside classifierBayes
  trainingSetInfo = models{nClasses+2};
end

% add extra things
trainingSetInfo.examplesWithLabel = examplesWithLabel;

% restore to make sure
models{nClasses+2} = trainingSetInfo; % kept for compatibility


%
% Create the "trainedClassifier" structure to be returned
%

% Eventually, a lot of the information inside "models" should be
% moved here, but we will keep the structure of "models" as is for now

% what the classifier code returns
trainedClassifier.models          = models;

% other information
trainedClassifier.trainingSetInfo = trainingSetInfo;
trainedClassifier.when            = datestr(now);
trainedClassifier.classifier      = classifierType;
trainedClassifier.classifierParameters      = classifierParameters;
trainingSetInfo;
