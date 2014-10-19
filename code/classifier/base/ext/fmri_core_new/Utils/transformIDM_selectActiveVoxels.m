% transformIDM_selectActiveVoxels(info,data,meta,nToKeep)
% Returns an IDM containing only the <nToKeep> most active voxels from each
% condition.  Thus, it keeps at most (#ROIS * #conditions * #nToKeep) voxels.
% Example:  
%  [info,data,meta] = transformIDM_selectActiveVoxels(info,data,meta,20);
%
% Dependencies:
% - transformIDM_selectActiveVoxels depends on mri_infoTrials
%
% History
% - 25 Aug 02 - fp - created from process_dataToExamples
%
%

function [info,ndata,nmeta] = transformIDM_selectActiveVoxels(varargin)
  info = varargin{1};
  data = varargin{2};
  meta = varargin{3};
  nToKeep = varargin{4};
  MergeOtherCond = 0;
  if (length(varargin) > 4)
    MergeOtherCond = varargin{5};
  end
  
% gather information about number,type and length of trials and other dataset parameters
% [ntrials,nvoxels,nconds,minTrialLenCond,ntrialsCond,trialBegin,trialEnd] = IDMinformation( info,data,meta,meta.study );
  IDM_information=IDMinformation( info,data,meta,meta.study );
  ntrials=IDM_information.nTrials;
  nvoxels=IDM_information.nVoxels;
  nconds=IDM_information.nConds;
  minTrialLenCond=IDM_information.minTrialLenCond;
  ntrialsCond=IDM_information.nTrialsCond;

  information = localInformation(meta.study);
  dimensions  = information.dimensions;
  dimx        = dimensions(1); dimy = dimensions(2); dimz = dimensions(3);
  
  % finds the top <nToKeep> most active voxels in each condition
  % and keeps all of those - at most #ROIS*#conditions*#nToKeep voxels
  % will be kept
  nvoxels = size(data{1},2); % WARNING: fails if data{1}==[]
  % a count of how often each voxel ends up in the top <nToKeep>,
  % across conditions
  activeCount = zeros(1,nvoxels);
  
  % get the Tvalues
  if MergeOtherCond
    [pvalueMaps] = mri_computeTvalues (info,data,meta,0);
  else
    [pvalueMaps] = mri_computeTvalues (info,data,meta);
  end

  nconds    = size(pvalueMaps,1);
  % calculate per condition ranks of voxels    
  for c=2:nconds
    pvalues   = pvalueMaps{1,c};
    % sort them by pvalue (smaller first) and get the top nToKeep
    ensemble  = sortrows(([pvalues;(1:nvoxels)])',1);
    topVoxels = ensemble(1:nToKeep,2);
    activeCount(topVoxels) = activeCount(topVoxels) + 1;
  end
  
  activeVoxels = find(activeCount > 0); % these will be used
  nactive      = length(activeVoxels);
  
  % go through data, crop to these voxels
  % and update the meta information to reflect it
  nmeta = meta;
  nmeta.nvoxels = nactive;
  nmeta.colToCoord = meta.colToCoord(activeVoxels,:);
  nmeta.coordToCol = zeros(size(meta.coordToCol));
  % could be made more efficient, but few voxels anyway
  % and this is clearer
  for v=1:nactive
    vc = nmeta.colToCoord(v,:);
    vx = vc(1); vy = vc(2); vz = vc(3);      
    nmeta.coordToCol(vx,vy,vz) = v;
  end
  
  ndata = cell(meta.ntrials,1);    
  for nt = 1:1:meta.ntrials
    ndata{nt} = data{nt}(:,activeVoxels);
  end
 
  clear data, meta;

    
