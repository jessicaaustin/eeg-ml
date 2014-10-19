% Builds a corresponding IDM out of an array of examples
% WARNING: works only for Naive Bayes and the categories dataset
% 
% This is used to invert the processing done by any of the idmToExample_<experiment>
% functions. These transform each trial in an IDM into an example.
% Each feature corresponds to a given voxel/"time point inside the trial"
% combination.
%
% You might want to do this if you have some vector of weights over
% features and want to see what it looks like on 3D, for instance.
% In that case, the examples array passed in would be a single row vector.
%
%
% In:
% - an example array (dimensions #examples * #features)
% - either
%   - an expInfo structure (if you are inverting an idmToExample)
%     (contains a Meta and a map of features to voxels+time within trial)
%   or
%   - a Meta + a featureToIDM map (not yet implemented) 
%
% Out:
% - an IDM
%
% Dependencies:
% - exampleToIDM depends on 
%
% History
% - 12 Sep 2002  - fp - created (it was about time!).


function [info,data,meta] = exampleToIDM( varargin )
  
  l = length(varargin);
  if l < 1
    fprintf(1,'syntax: exampleToIDM(examples,expInfo)\n\n');
    return;
  end
 
  examples     = varargin{1};
  expInfo      = varargin{2};
  meta         = expInfo.meta;
  experiment   = expInfo.experiment;
  featureToIDM = expInfo.featureToIDM;
  %  dataType   = expInfo.dataType; % not there in a few of the past saves

  % figure out more about the IDM
  featureVoxels = featureToIDM(:,1);
  featureTimes  = featureToIDM(:,2); 
  voxels        = sort( unique( featureVoxels ));
  timePoints    = sort( unique( featureTimes ));
  minTrialLen   = length(timePoints);
  nVoxels       = length(voxels);
  [nTrials,nFeatures] = size(examples);
  
  % now create the IDM structures
  for nt=1:nTrials
    info(nt).len  = minTrialLen;
    info(nt).cond = -1;
    data{nt} = zeros(minTrialLen,nVoxels);
  end
  
  % precompute a few useful things
  % assumes all voxels appear at all time steps
  for t=1:minTrialLen
    % find features with this timestamp
    featureNumbers{t} = find(featureTimes==t);
    % find the corresponding voxels
    voxelNumbers{t}   = featureVoxels(featureNumbers{t});
    ensemble          = sortrows([featureNumbers{t},voxelNumbers{t}],2);  
    featureOrder{t}   = ensemble(:,1);
  end
  
  % and convert the features into data
  for nt=1:nTrials
    featureVector  = examples(nt,:);
    for t=1:minTrialLen
      data{nt}(t,:)  = featureVector(featureOrder{t});
    end
  end

  
  
  
  
function [] = testThis()
  
% create a dummy dataset
  data{1} = randn(3,10);
  data{2} = randn(4,10);
  data{3} = randn(3,10);
  data{4} = randn(2,10);
  
  info(1).len = 3;
  info(1).cond = 3;
  info(2).len = 4;
  info(2).cond = 2;
  info(3).len = 3;
  info(3).cond = 3;
  info(4).len = 2;
  info(4).cond = 1;
  
  meta.study = 'categories';
  dataType = 'unfolded';
  
  % create examples
  [examples,labels,expInfo] = idmToExamples_1ofN(info,data,meta,dataType);
  
  % now rebuild IDM from them
  [ninfo,ndata,nmeta] = exampleToIDM( examples,expInfo );  
  
  % and check
  for i=1:3
    isequal(data{i}(1:11,:),ndata{i})
  end