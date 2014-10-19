% transformIDM_joinTrials(info,data,meta) 
%
% Transforms info, data, and meta by concatenating all trials into a single
% trial.  The returned info has info.cond=-1 to indicate there is no single
% condition for the merged trial.  It also contains info.snapshotCond(i) set
% to the condition of the trial from which the i-th snapshot was
% taken.
%
% Access the single trial data via data{1}, and the info for this trial via
% info(1).
%
% Example: 
% [info2,data2,meta2] = transformIDM_mergeTrials(info,data,meta) 
%
% History
% - 8/21/02 TMM Created file.
% - 27 Nov 02 - fp - reworked it to make it faster

function [rinfo,rdata,rmeta] = transformIDM_joinTrials(info,data,meta)

  ntrials = length(data);
  nvoxels = size(data{1},2);
 
  trialSequence = [info(:).cond];
  trialLens     = [info(:).len];
  totalLen      = sum(trialLens);
 
  rdata = cell(1,1);
  rdata{1} = zeros(totalLen,nvoxels);
  rinfo(1).len  = totalLen;
  rinfo(1).cond = -1;
  rmeta = meta;
  rmeta.ntrials = 1;
  
  idx = 1;

  for nt=1:ntrials
    len = info(nt).len;
    rdata{1}(idx:idx+len-1,:) = data{nt};    
    idx = idx + len;
  end;
  
