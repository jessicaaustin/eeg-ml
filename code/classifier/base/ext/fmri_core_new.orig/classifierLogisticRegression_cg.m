% Train a logistic classifier
% 
% In:
% - trainSet       (#trainexamples*#features) matrix
% - trainLabels    (#trainexamples)   column vector
% - classifierParameters it should be a cell, which
% contains two elements: first- stopCriterion and second- lamda.
% If this cell is empty, the defaults are:
%      stopCriterion = 0.001;
%      lamda         = 10;
%
% Out:
% - models - contains learnt models
%
% Dep:
%
% History: 
% created May 02, 2005 by Wei Wang. 
% optimization was adapted from Francisco Pereira.
%
% Known bugs:
%
% Ex:
% - [models]=classifierLogisticRegression_cg([1.0 1.5; 1.1 2.0; 2.0 1.0; 2.1 2.0],[0;0;1;1])
%
% Reference: 
% - Machine learning by Tom Mitchell
% - optimization method uses conjugate gradient
% - part of optimization codes was adapted from FP
%  

function [models] = classifierLogisticRegression_cg( varargin )
  
 l = length(varargin);
  if l < 3
    fprintf('syntax: classifierLogisticRegression_cg(trainSet,trainLabels, parameters)\n');
    return;
  elseif l > 3
    fprintf('syntax: classifierLogisticRegression_cg(trainSet, trainLabels, parameters)');
    return;
  end
  
  
  trainSet    = varargin{1};
  trainSet    = [ones(size(trainSet,1),1) trainSet];
  trainLabels = varargin{2};
  classifierParameters = varargin{3};
  if length(classifierParameters) > 2
    fprintf('syntax: parameters for classifierLogistiRegression should be 2');
    return;
  end
  
  nTrain      = size(trainSet,1);
  nFeatures   = size(trainSet,2);
  nLabels     = size(trainLabels,2);
  
  if length(classifierParameters) == 0
      stopCriterion = 0.001;
      lamda         = 10;
  else
      stopCriterion = classifierParameters{1};
      lamda         = classifierParameters{2};
  end

  
  if nTrain == 0
    % not very graceful, but meta experiment dummy uses this
    models = {}; return;
  end
    
  sortedLabelValues = sort(unique(trainLabels));
  nClasses          = length(sortedLabelValues);

  % each column is the weight for one class
  weights = zeros(nFeatures,nClasses-1);
  
  % estimate weights
  preLogL = -Inf;
  logistic_g = zeros(nFeatures,nClasses-1);
  
  % first round
  [grad, py_k, delta] = gradient(sortedLabelValues,trainSet,trainLabels,weights);
  logistic_g = grad - lamda*weights;
      
  logistic_g=reshape(logistic_g,nFeatures*(nClasses-1),1);
  logistic_d=-logistic_g;
      
  logistic_Q=logisticComputeQ(trainSet, py_k)-lamda*eye(nFeatures*(nClasses-1));
  logistic_alpha=-1*logistic_g'*logistic_g/(logistic_g'*logistic_Q*logistic_d);
  weights = weights-logistic_alpha*reshape(logistic_d,nFeatures,nClasses-1);
  % n-1 rounds
  logistic_gg=zeros(nFeatures,nClasses-1);
  
  %for nDirections=1:(nFeatures*(nClasses-1)-1)
  nIterations = 0;
  while 1
          [grad, py_k, delta] = gradient(sortedLabelValues,trainSet,trainLabels,weights);
          logL = compute_logLikelihood(py_k,delta);
          if logL > preLogL + stopCriterion
              preLog = logL;
          else
              break;
          end
          logistic_gg = grad - lamda*weights;
          logistic_gg=reshape(logistic_gg,nFeatures*(nClasses-1),1);
          logistic_beta=logistic_gg'*logistic_gg/(logistic_g'*logistic_g);
          logistic_d=-logistic_gg+logistic_beta*logistic_d;
          logistic_Q=logisticComputeQ(trainSet, py_k)-lamda*eye(nFeatures*(nClasses-1));
          logistic_alpha=-1*logistic_gg'*logistic_gg/(logistic_gg'*logistic_Q*logistic_d);
          weights = weights-logistic_alpha*reshape(logistic_d,nFeatures,nClasses-1);
          logistic_g=logistic_gg;
          nIterations = nIterations + 1;
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
  %trainingSetInfo.classPriors      = classPriors;
  models{nClasses+2} = trainingSetInfo;
  models{nClasses+3} =[];
  
function [grad, py_k, delta]=gradient(sortedLabelValues,trainSet,trainLabels,weights)
    nClasses=length(sortedLabelValues);
    [nTrain,nFeatures] = size(trainSet);
    
    weights=[weights zeros(size(weights,1),1)];  
    tmp=exp(trainSet*weights);
    dsum=sum(tmp,2);
  
    py_k=tmp ./ repmat(dsum,1,nClasses);
    delta=repmat(trainLabels,1,nClasses)==repmat(sortedLabelValues',nTrain,1);
      
    stepk=zeros(size(weights));
    for k=1:(nClasses-1)
          errork=delta(:,k)-py_k(:,k);
          stepk(:,k)=sum(trainSet .* repmat(errork,1,nFeatures), 1)';
    end
    grad = stepk(:,1:(end-1)); 
    py_k = py_k(:,1:(end-1));
      
    
function logistic_Hessian=logisticComputeQ(trainSet,py_k)
    [nTrain, nFeatures]=size(trainSet);
    [nTrain, nClassesLessOne]=size(py_k);
    
    logistic_Hessian=zeros(nClassesLessOne*nFeatures);
    
    for i = 1:nTrain
        x=trainSet(i,:)';
        tmp_X = repmat(x*x',nClassesLessOne, nClassesLessOne);
        p=py_k(i,:);
        p=repmat(p,nFeatures,1);
        p=reshape(p,nClassesLessOne*nFeatures,1);
        tmp_P = p*p';
        logistic_Hessian = logistic_Hessian + tmp_X.*tmp_P;
    end
    
function [L] = compute_logLikelihood(py_k,delta)
    py_k_last = ones(size(py_k,1),1) - sum(py_k,2);
    py_k = [py_k py_k_last];
    L_examples = log( sum(delta .* py_k,2) );
    L = sum(L_examples,1);

