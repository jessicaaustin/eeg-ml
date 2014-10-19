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
% - [models]=classifierSMLR([1.0 2.0; 1.1 2.0; 2.0 1.0; 2.1 2.0],[1;1;0;0],{})
%
% Reference: 
% - Machine learning by Tom Mitchell
%   This funciton uses the steepest descent with the automated stopping
%   rule as the optimization method.  

function [models] = classifierSMLR( varargin )
  
  l = length(varargin);
  if l < 3
    fprintf('syntax: classifierSMLR(trainSet,trainLabels, parameters)\n');
    return;
  elseif l > 3
    fprintf('syntax: classifierSMLR(trainSet, trainLabels, parameters)');
    return;
  end
  
  
  trainSet    = varargin{1};
  trainLabels = varargin{2};
  classifierParameters = varargin{3};
  
  [trainSet, mean_trainLabels] = correlation_train(trainSet, trainLabels);
  [models] = classifierSMLR_core(trainSet, trainLabels, classifierParameters);
  [models] = changeWeight(models, mean_trainLabels);
  
  
function [new_d, mean_d] = correlation_train(data,labels)
    conds=sort(unique(labels));
    numConds=length(conds);
    [numRows, numCols]=size(data);
    mean_d=zeros(numConds,numCols);
    for j=1:numRows
        k=find(conds==labels(j));
        mean_d(k,:)=mean_d(k,:)+data(j,:)/sum(labels==labels(j));
    end
    
    new_d = zeros(numRows, numCols);
    for j=1:numRows
        k=find(conds==labels(j));
        new_d(j,:) = (data(j,:)-mean(data(j,:))).*(mean_d(k,:)-mean(mean_d(k,:)))/(std(data(j,:))*std(mean_d(k,:)));
    end
        
function [models_new] = changeWeight(models, mean_d)
    
    lastpos           = size(models,1);
    weights=models{lastpos-2};
    [nFeatures, nClasses]=size(weights);
    
    mean_d=[ones(size(mean_d,1),1) mean_d];
    c=mean_d';
    
    c=c(:,1:nClasses);
    weights_new=zeros(nFeatures,nClasses);
    for j=1:nClasses
        weights_new(:,j)=(c(:,j)-mean(c(:,j))).*weights(:,j)/std(c(:,j));
    end
    models_new = models;
    models_new{lastpos-2}=weights_new;
    models_new{lastpos-1}.uncombinedWeights=weights;


