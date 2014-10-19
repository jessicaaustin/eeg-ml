% transformIDM_selectActiveVoxels(info,data,meta,nToKeep)
%
% Returns an IDM containing only the <nToKeep> most active voxels
% across selected conditions (default: all but fixation).
%
% Uses summarizeELE_rankByMeanTest to rank the voxels, and 
% the procedure is as follows:
%
% - Find p-values for a t-test of every voxel against 0,
% separately for each condition, and rank by  mean != 0 or mean > 0
%
% - Use the ranks to select voxels
% Conceptually, this works by taking the per condition rankings
% Cond  2   ... last condition being considered
%       f(1)    f(1)
%       f(2)    f(2)
% and  -------------- sliding a window down until the number left
% above is the # to keep.
%
% Note that voxels near the top of several of the rankings will
% only get selected once.
% 
% In:
% - info+data+meta
% - # of voxels to keep
% - optionals
%   - conditions - conditions to use for ranking (default is all)
%   - test type - either 'greaterThanZero' (default) or 'differentFromZero'
%
% Out:
% - an IDM containing only the selected voxels

% Dependencies:
%
% History
% - 10 Nov 2004 - fpereira - created from prior code
%

function [info,data,meta,voxelsToKeep] = transformIDM_selectActiveVoxV2(varargin)
  
  l = length(varargin);
  
  if l < 4; help transformIDM_selectActiveVoxactV2; return
  else    
    info    = varargin{1};
    data    = varargin{2};
    meta    = varargin{3};
    nToKeep    = varargin{4};

    trialSequence = [info(:).cond];
    trialLens     = [info(:).len];
    taskTrials    = find(trialSequence > 1);
    conditions    = unique(trialSequence(taskTrials));
    testType      = 'greaterThanZero';
    if l > 4
      conditions = varargin{5};
      if l > 5
	testType = varargin{6};
      end
    end
    nConds        = length(conditions);
    nVoxels       = size(data{1},2);
    nTrials       = length(trialSequence);
  end

  if nToKeep > nVoxels
    if nVoxels
      fprintf('warning: roi %s has %d voxels and you asked me to keep %d.\n',meta.roi,nVoxels,nToKeep);
      nToKeep = nVoxels;
    else
      fprintf('ERROR: roi %s does not contain voxels\n',meta.roi);
      return;
    end    
  end
  
  fprintf('transformIDM_selectActiveVoxact: keeping %d voxels from condition(s) ',nToKeep);
  fprintf('%d ',conditions); fprintf('\n');

  %% create a large matrix with all the condition trials, pass that
  taskTrialLens = trialLens(taskTrials);
  taskMatrix    = zeros(sum(taskTrialLens),nVoxels);
  labelMatrix   = zeros(sum(taskTrialLens),1);

  idx = 1;
  for nt = 1:nTrials
    if info(nt).cond > 1
      taskMatrix(idx:idx+info(nt).len-1,:) = data{nt};
      labelMatrix(idx:idx+info(nt).len-1)  = info(nt).cond;
      idx = idx + info(nt).len;
    end
  end
  
  [sortedVoxels] = summarizeELE_rankByMeanTest(taskMatrix,labelMatrix,testType);  
  voxelsToKeep   = sortedVoxels(1:nToKeep);
  
  %% select voxels
  [info,data,meta] = transformIDM_selectVoxelSubset(info,data,meta,voxelsToKeep);

  return;

  
function [] = testThis()
  
  subject = '05868';
  study   = 'data-categories';
  experiment = '1ofN';

  information = mri_reference(study);
  rois        = information.rois;
  
  % first load data
  rois = {'LIT'};
  %  [info,data,meta] = loadSubjectdata(study,subject,rois);
  [info,data,meta] = loadSubjectdata(study,subject,rois);
  % Turn each trial into a separate condition
  % (note: don't do this with 05675/05695, as they have missing trials)
  [info,data,meta] = transformIDM_separateBlocks(info,data,meta);
  % Transform the IDMs again, by unfolding them.
  [info,data,meta] = transformIDM_unfold(info,data,meta);

  % get information about IDM
  [nTrials,nVoxels,nConds,minTrialLenCond,ntrialsCond] = mri_infoTrials( info,data,meta,study);

  transformIDM_selectActiveVoxact(info,data,meta,20);
  
  
  % mock example set
  % - first 5 features  (nonactive->active) in condition 1
  % - second 5 features (nonactive->active) in condition 2
  part1 = randn(20,5);
  part1 = part1 + ones(20,1) * [0 0.25 0.5 0.75 1];
  part2 = randn(20,5);
  part2 = part2 + ones(20,1) * [0 0.25 0.5 0.75 1];
  examples = randn(40,10);
  examples(1:20,1:5)  = part1;
  examples(21:40,6:10) = part2;
  labels   = zeros(40,1);
  labels(1:20) = 1; labels(21:40) = 2;

  % make an IDM
  info(1).cond = 2; info(1).len = 20;
  data{1} = examples(1:20,:);
  info(2).cond = 3; info(2).len = 20;
  data{2} = examples(21:40,:);
  meta.dimx = 20; meta.dimy = 20; meta.dimz = 20;
  meta.colToCoord(1:10,:) = [(1:10)',ones(10,1),ones(10,1)];
  c = meta.colToCoord;
  indices = sub2ind([20,20,20],c(:,1),c(:,2),c(:,3));
  meta.coordToCol = zeros(20,20,20);
  meta.coordToCol(indices) = 1:10;
  meta.nvoxels = 10;
