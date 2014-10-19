% Average all the trials in each condition, yielding an IDM with
% one trial per condition.
%
% The average trial for each condition will be as long as the shortest
% of the trials of that condition that went into the average.
%
% In:
% - info, data and meta
%
% Out:
% - avginfo+meta (info and meta for the average trial data
% structure)
% - dataMean - data structure with the average trials
% - dataStdv - data structure with the standard deviation of each
% average computation
% 
% Notes:
% - condition 0 is ignored
% - any info fields other than cond and len are destroyed
%
% Examples:
% [avginfo,dataMean,meta,dataStdv] = transformIDM_avgTrialCondition(info,data,meta);
%
% History:
% - 17 Jul 2005 - fpereira - modified from previous code

function [avginfo,dataMean,meta,dataStdv] = transformIDM_avgTrialCondition( varargin )
  
%% process parameters

l = length(varargin);
if l < 1; fprintf('syntax: transformIDM_dataMeanByCond2(info,data,meta)\n');return;end

info = varargin{1}; data = varargin{2}; meta = varargin{3};
    
trialSequence = [info(:).cond];
trialLens     = [info(:).len];
nonzeroTrials = find(trialSequence > 0);
nvoxels       = size(data{1},2);
  
conditionsPresent = unique(trialSequence(nonzeroTrials));
nConditions       = length(conditionsPresent);
  
%% compute the average trial

dataMean = cell(nConditions,1);
dataStdv = cell(nConditions,1);  
  
for c = 1:nConditions
    
  % find the trials with this cond
  cond = conditionsPresent(c);  
  trialsOfThisCond  = find(trialSequence == cond);
  lensOfThoseTrials = trialLens(trialsOfThisCond);
  minlen  = min(lensOfThoseTrials);
  ntrials = length(trialsOfThisCond);
    
  % copy them into a 3D matrix, cropping them to the minimum trial length
  tmp = zeros(minlen,nvoxels,ntrials);
    
  for nt = 1:ntrials
    idx = trialsOfThisCond(nt);
    tmp(:,:,nt) = data{idx}(1:minlen,:);
  end

  % compute average and standard deviation trial for this condition
  dataStdv{c} = std(tmp,0,3);
  dataMean{c} = mean(tmp,3);

  % fill in the info structure 
  avginfo(c).cond = cond;
  avginfo(c).len  = minlen;
end
  
% finish up the meta information
meta.ntrials = nConditions;
  

