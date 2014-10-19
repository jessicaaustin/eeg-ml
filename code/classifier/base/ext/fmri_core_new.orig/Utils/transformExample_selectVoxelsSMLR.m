% transformExample_selectVoxelsSMLR(examples,nToKeep,classifier,selMeth)
%
% Returns an examples containing only the <nToKeep> most active voxels
% by the SMLR weights.
%
% selMethod:
%   1.  mean of abs of weights over all possible sets of weights
%   2.  max of abs of weights over all possible sets of weights
%   3.  diff of weights over all possible sets of weights
% Example:  
%  [rexamples,rclassifier] = transformExample_selectVoxelsSMLR(examples,20,classifier,1);
%
% Dependencies:
% - 
%
% History
% - 25 May 05 - Wei - created 

function [rexamples, rclassifier, idxVox] = transformExample_selectVoxelsSMLR(varargin)
  
  l = length(varargin);
  
  if l < 4
    fprintf('syntax: transformExample_selectVoxelsSMLR(examples,<# to keep>,classifier,method)\n');
    return;
  end
    
  examples = varargin{1};
  nToKeep = varargin{2};
  classifier = varargin{3};
  selection = varargin{4};
    
  models = classifier.models;
  lastpos = size(models,1);
  % use uncombinded weights to pick the voxels
  weights_all=models{lastpos-1}.uncombinedWeights;
  weights=weights_all(2:end,:);
  
  fprintf('transformExample_selectVoxelsSMLR: keeping %d voxels\n',nToKeep);
  
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
  [tmp, idx] = sort(tmp,'descend');
  activeVoxes =idx(1:nToKeep);
  
  % output the combined weights
  weights=models{lastpos-2};
  rclassifier = classifier;
 
  rexamples = examples(:,activeVoxes);
  rweights = [weights(1,:); weights(activeVoxes,:)];
  rclassifier.models{lastpos-2} = rweights;
  idxVox = idx;