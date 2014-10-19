% transformExample_selectVoxelsLR(examples,nToKeep,classifier,selMeth)
%
% Returns an examples containing only the <nToKeep> most active voxels
% by the logisticRegression weights.
%
% selMethod:
%   1.  mean of abs of weights over all possible sets of weights
%   2.  max of abs of weights over all possible sets of weights
%   3.  diff of weights over all possible sets of weights
%
% Input:
%   examples -- a matrix with dimension of (# of instances)*(# of features)
%   nToKeep  -- a integer greater than 0
%   classifier -- a learnt classifie
%   setMeth -- 1/2/3, see above
% Output:
%   rexamples -- a matrix with of (# of instances)*(# of selected
%   features), which contains the only seleted features
%   rclassifer -- same as input except the weight has been selected,
%   selected weights are associated with the seleted features
%   idxVox -- a (# of features)*1 array which indicates the index of features
%   according to the importance in the descend order
%
% Example:  
%  [rexamples,rclassifier, idx] = transformExample_selectVoxelsLR(examples,20,classifier,1);
%
% Dependencies:
% - 
%
% History
% - 06 Jul 05 - Wei - created 

function [rexamples, rclassifier, tmp] = transformExample_selectVoxelsLR(varargin)
%function [rexamples, rclassifier, idxVox] = transformExample_selectVoxelsLR(varargin)
  
  l = length(varargin);
  
  if l < 4
    fprintf('syntax: transformExample_selectVoxelsLR(examples,<# to keep>,classifier,method)\n');
    return;
  end
    
  examples = varargin{1};
  nToKeep = varargin{2};
  classifier = varargin{3};
  selection = varargin{4};
    
  models = classifier.models;
  lastpos = size(models,1);

  weights_all=models{lastpos-2};
  weights=weights_all(2:end,:);
  
  fprintf('transformExample_selectVoxelsLR: keeping %d voxels\n',nToKeep);
  
  switch selection
      % average of weights of each voxels
      case 1
          tmp = mean(abs(weights),2);
          
      % max over the classes
      case 2
          tmp = max(abs(weights),[],2);
    
      % max-min over classes
      case 3
          tmp = max(weights,[],2) - min(weights,[],2);
  end
%  keyboard
  %[tmp, idx] = sort(tmp,'descend');
  %[tmp,idx]=sort(tmp);
  idx=tmp;
  activeVoxes =idx(1:nToKeep);
  
  rclassifier = classifier;
 rexamples=[]
 % rexamples = examples(:,activeVoxes);
 % rweights = [weights_all(1,:); weights(activeVoxes,:)];
 % rclassifier.models{lastpos-2} = rweights;
  idxVox = idx;