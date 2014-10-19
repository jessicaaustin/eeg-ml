% Given an IDM, creates examples for a fixation learning task
%
% Essentially, turns a trial (#voxels over a time interval
% corresponding to performance of the task in a given condition)
% into an example.
%
% In:
% - IDM
%
% Out:
% - examples, labels and experimental information
%
% Dependencies:
% - expCategories_fixation depends on mri_infoTrials
%
% History:
% - 09 Aug 02 - fp - create from exp_categories_fixation in fmri_devel
% - 28 Aug 02 - tm - removed fourth argument (dataType) because it was never used!
%
%
  
function [examples,labels,expInfo] = idmToExamples_fixation( varargin )
  
  l = length(varargin);
  if l < 3
    fprintf(1,'syntax: idmToExamples_fixation(info,data,meta)');
    examples = []; labels = [];
    return;
  else
    info = varargin{1};
    data = varargin{2};
    meta = varargin{3};
    if l > 3; dataType = varargin{4}; else dataType = 'dummy'; end
    study = meta.study;
  end
  
  % gather information about number,type and length of trials
  %[ntrials,nvoxels,nconds,minTrialLenCond,ntrialsCond] = mri_infoTrials( info,data,meta,study );
  IDM_information=IDMinformation( info,data,meta,meta.study );
  ntrials=IDM_information.nTrials;
  nvoxels=IDM_information.nVoxels;
  nconds=IDM_information.nConds;
  minTrialLenCond=IDM_information.minTrialLenCond;
  ntrialsCond=IDM_information.nTrialsCond;
  trialSequence=IDM_information.trialSequence;
  
  minTrialLen = min(minTrialLenCond(1:1:nconds));
  maxCond     = max(trialSequence);
    
  %
  % set up for outputting ex  %
  % Examples are stored as 
  % v1...vn (time1) ... v1 ... vn (time t)
  %  
  nfeatures = minTrialLen * nvoxels;
  nexamples = sum(ntrialsCond(:))-ntrialsCond(1);  

  examples  = zeros(nexamples,nfeatures);
  labels    = zeros(nexamples,1);

  % This will keep information to allow mapping from feature space
  % to IDM space, i.e. given a feature we can say which voxel,trial
  % and time it corresponds to. This is useful for taking any
  % vector of values in feature space (such as weights from a
  % learnt classifier and seeing which voxels and times they
  % correspond to.
  % expInfo.featureToIDM is a matrix with dimensions #features * 2
  % The row for each feature contains [voxel column, time within trial]  

  expInfo.experiment   = 'fixation';
  expInfo.dataType     = dataType;
  expInfo.meta         = meta;
  expInfo.featureToIDM = zeros(nfeatures,2);
  
  % Examples are stored as v1...vn (time1) ... v1 ... vn (time t)
  % so reproduce this structure for this case
  
  ftoidm = zeros(nvoxels,2);
  ftoidm(:,1) = (1:nvoxels)';  
  featureIdx = 1;
  
  for t = 1:1:minTrialLen
    featureIdxNext = featureIdx + nvoxels;	
%   examples( exampleIdx, featureIdx:(featureIdxNext-1) ) = data{nt}(t,:);
  
    ftoidm(:,2) = ones(nvoxels,1)*t;
    expInfo.featureToIDM( featureIdx:(featureIdxNext-1),: ) = ftoidm;
    
    featureIdx = featureIdxNext;
  end

  % 
  % Finally, create the examples
  %
  
  exampleIdx = 1;
  for nt = 1:1:ntrials
    len  = info(nt).len;
    cond = info(nt).cond;

    if cond > 0 & cond <= maxCond
      featureIdx = 1;
  
      % concatenate the signal at all voxels for each time point
      for t = 1:1:minTrialLen
	featureIdxNext = featureIdx + nvoxels;	
	examples( exampleIdx, featureIdx:(featureIdxNext-1) ) = data{nt}(t,:);

	featureIdx = featureIdxNext;
      end

      % Decide what to do with this example based on condition #
      % The decision could be based on anything else, such as
      % the stimulus, a behavioural indicator, etc. All this information
      % is available from info.
      
      if cond > 1
	labels(exampleIdx,1) = 2;
      else
	labels(exampleIdx,1) = 1;
      end
      
      exampleIdx = exampleIdx + 1;
    else
      % just ignore condition 0
    end
  end

  
