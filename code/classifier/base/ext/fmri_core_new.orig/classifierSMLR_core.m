% Train a SMLR classifier
% 
% In:
% - trainSet       (#trainexamples*#features) matrix
% - trainLabels    (#trainexamples)   column vector
% - classifierParameters it should be a cell, which
% contains three elements: first- initial stepsize and second - stopCriterion
% and third - lamda. If this cell is empty, the defaults are:
%      stepsize      = 0.01;
%      stopCriterion = 0.001;
%      lamda         = 10;
%
% Out:
% - models - contains learnt models
%
% Dep:
%
% History: 
% created Mar 20, 2005 by Wei Wang. 
% optimization was adapted from Francisco Pereira.
%
% Known bugs:
%
% Ex:
% - [models]=classifierSMLR_core([1.0 2.0; 1.1 2.0; 2.0 1.0; 2.1 2.0],[1;1;0;0],{})
%
% Reference: 
% - Machine learning by Tom Mitchell
%  

function [models] = classifierSMLR_core( varargin )
  
  l = length(varargin);
  if l < 3
    fprintf('syntax: classifierSMLR_core(trainSet,trainLabels, parameters)\n');
    return;
  elseif l > 3
    fprintf('syntax: classifierSMLR_core(trainSet, trainLabels, parameters)');
    return;
  end
  
  trainSet    = varargin{1};
  trainSet    = [ones(size(trainSet,1),1) trainSet];
  trainLabels = varargin{2};
  classifierParameters = varargin{3};
  if length(classifierParameters) > 3
    fprintf('syntax: parameters for classifierSMLR_core should be 3');
    return;
  end
  
  nTrain      = size(trainSet,1);
  nFeatures   = size(trainSet,2);
  nLabels     = size(trainLabels,2);
  
  if length(classifierParameters) == 0
      stepsize      = 0.01;
      stopCriterion = 0.001;
      lamda         = 10;
  else
      stepsize      = classifierParameters{1};
      stopCriterion = classifierParameters{2};
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
  
  nIterations=0;
  preLogL = -Inf;
  
  while 1
      
      tmp=exp(trainSet*weights);
      dsum=sum(tmp,2);
  
      py_k=tmp ./ repmat(dsum,1,nClasses);
      delta=repmat(trainLabels,1,nClasses)==repmat(sortedLabelValues',nTrain,1);
      
      stepk=zeros(size(weights));
      for k=1:(nClasses)
          errork=delta(:,k)-py_k(:,k);
          stepk(:,k)=sum(trainSet .* repmat(errork,1,nFeatures), 1)';
      end
      
      weights=weights + stepsize*(stepk - lamda*weights);
      
      logL = compute_logLikelihood(py_k,delta);
      if ~isnan(logL)
          if logL > preLogL + stopCriterion
              preLogL = logL;
          elseif stepsize > 0.00001
              stepsize = stepsize/2;
              preLogL = logL;
          else
              break;
          end
      else
          stepsize = stepsize / 10;
          weights = zeros(nFeatures,nClasses);
          nIterations = 0;
          fprintf('decimating stepsize as %d and restarting\n', stepsize);
      end
      
      nIterations = nIterations + 1;
      %fprintf(' nIterations %d\n LogL %d\n stepsize %d\n',nIterations,logL,stepsize);
  end
  fprintf('running %d iterations\n', nIterations);
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
  trainingSetInfo.classPriors       = zeros(nClasses,1);
  for l=1:nClasses
    trainingSetInfo.classPriors(l) = length(find(trainLabels==sortedLabelValues(l)));
  end
  trainingSetInfo.classPriors = trainingSetInfo.classPriors/nTrain;   
  
  models{nClasses+2} = trainingSetInfo;
  models{nClasses+3} =[];
  
  
function [L] = compute_logLikelihood(py_k,delta)

    L_examples = log( sum(delta .* py_k,2) );
    L = sum(L_examples,1);


