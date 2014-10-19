% Train a logistic classifier
% 
% In:
% - trainSet       (#trainexamples*#features) matrix
% - trainLabels    (#trainexamples)   column vector
% - classifierParameters it should be a cell, which
% contains two elements: first-  stepsize and second - nIterations
% and third - lamda. If this cell is empty, the defaults are:
%      stepsize      = 0.00001;
%      nIterations   = 10000;
%      lamda         = 10;
%
% Out:
% - models - contains learnt models
%
% Dep:
%
% History: 
% created Mar 02, 2005 by Wei Wang. 
%
% Known bugs:
%
% Ex:
% - [models]=classifierLogisticRegression_constant([1.0 2.0; 1.1 2.0; 2.0 1.0; 2.1 2.0],[1;1;0;0], {})
%
% Reference: 
% - Machine learning by Tom Mitchell
%   This funciton uses the steepest descent with the constant stepsize as the optimization method.  

function [models] = classifierLogisticRegression_constant( varargin )
  
  l = length(varargin);
  if l < 3
    fprintf('syntax: classifierLogisticRegression_constant(trainSet,trainLabels, parameters)\n');
    return;
  elseif l > 3
    fprintf('syntax: classifierLogistiRegression_constant(trainSet, trainLabels, parameters)');
    return;
  end
  
  
  trainSet    = varargin{1};
  trainSet    = [ones(size(trainSet,1),1) trainSet];
  trainLabels = varargin{2};
  classifierParameters = varargin{3};
  if length(classifierParameters) > 3
    fprintf('syntax: parameters for classifierLogistiRegression_constant should be 3');
    return;
  end
  
  nTrain      = size(trainSet,1);
  nFeatures   = size(trainSet,2);
  nLabels     = size(trainLabels,2);
  
  if length(classifierParameters) == 0
      stepsize      = 0.00001;
      nIterations   = 10000;
      lamda         = 10;
  else
      stepsize      = classifierParameters{1};
      nIterations   = classifierParameters{2};
      lamda         = classifierParameters{3};
  end
  
  if nTrain == 0
    % not very graceful, but meta experiment dummy uses this
    models = {}; return;
  end
    
  sortedLabelValues = sort(unique(trainLabels));
  nClasses          = length(sortedLabelValues);

  % each column is the weight for one class
  weights = zeros(nFeatures,nClasses);
  
  % estimate weights
  for step=1:nIterations
      weights(:,nClasses)=zeros(nFeatures,1);
      
      tmp=exp(trainSet*weights);
      dsum=sum(tmp,2);
  
      py_k=tmp ./ repmat(dsum,1,nClasses);
      delta=repmat(trainLabels,1,nClasses)==repmat(sortedLabelValues',nTrain,1);
      
      stepk=zeros(size(weights));
      for k=1:(nClasses-1)
          errork=delta(:,k)-py_k(:,k);
          stepk(:,k)=sum(trainSet .* repmat(errork,1,nFeatures), 1)';
      end
      
      weights=weights + stepsize*(stepk - lamda*weights);
  end
  
  
  %% Now prepare output in a cell array
  models = cell(nClasses+1,1);

  % Generative model - each cell contains a cell array
  % where each cell has one parameter - mean, covariance matrix, etc

  % Discriminative model - a cell array of sets of weight
  models{nClasses+1} = weights;
  
  % Training Set information
  trainingSetInfo.nExamples         = nTrain;
  trainingSetInfo.nFeatures         = nFeatures;
  trainingSetInfo.nLabels           = nLabels;
  trainingSetInfo.nClasses          = nClasses;
  trainingSetInfo.sortedLabelValues = sortedLabelValues;
  %trainingSetInfo.classPriors      = classPriors;
  models{nClasses+2} = trainingSetInfo;
  models{nClasses+3} =[];
  
    

