% Gather information about number, type and length of trials
%
% In:
% - info+data+meta
%
% Out:
% - an "information" structure, with the following fields
% - ntrials - # trials
% - nvoxels - # voxels
% - nconds  - # conditions
% - minTrialLenCond   - smallest trial length in each condition >= 1
% - nTrialsCond       - # of trials of each condition >= 1
% - trialSequence     - sequence of trial conditions
% - trialLengths      - sequence of trial lengths
% - taskTrials        - indices of trials in "trialSequence" with condition >= 1
% - taskTrialSequence - conditions of those trials
% - taskTrialLengths  - lengths    of those trials

% Notes:
%
% - 23 Mar 2005 - fpereira - created

function [information] = IDMinformation( varargin )

%
% Process parameters
%

l = length(varargin);
if l < 3; help IDMinformation; return; end
info = varargin{1};
data = varargin{2};
meta = varargin{3};
if l > 4; study = varargin{4}; else; study = meta.study; end

ntrials = max(size(data));
nvoxels = size(data{1},2);

trialSequence = [info(:).cond];
trialLengths  = [info(:).len];
taskTrials    = find(trialSequence > 0);
taskTrialSequence = trialSequence(taskTrials);
taskTrialLengths  = trialLengths(taskTrials);

conditions  = unique( taskTrialSequence );
nconds      = length(conditions);

conditionToTrial = cell(nconds,1);

for c = 1:nconds
  conditionNumber = conditions(c);

  % trials in this condition
  conditionToTrial{c} = find(trialSequence == conditionNumber);

  % minimum length of a trial in this condition
  minTrialLenCond(c) = min( trialLengths( conditionToTrial{c}) );

  % # of trials in this condition
  nTrialsCond(c) = length(conditionToTrial{c});
end


% pack results

information.nTrials = ntrials;
information.nVoxels = nvoxels;
information.nConds  = nconds;
information.minTrialLenCond = minTrialLenCond;
information.nTrialsCond     = nTrialsCond;
information.trialSequence = trialSequence;
information.trialLengths  = trialLengths;
information.taskTrials    = taskTrials;
information.taskTrialSequence = taskTrialSequence;
information.taskTrialLengths  = taskTrialLengths;
information.conditionToTrial = conditionToTrial;
