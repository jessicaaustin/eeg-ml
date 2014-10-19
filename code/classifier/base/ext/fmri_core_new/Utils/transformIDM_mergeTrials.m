% transformIDM_mergeTrials(info,data,meta) 
%
% Transforms info, data, and meta by concatenating all trials into a single
% trial.  The returned info has info(1).cond=-1 to indicate there is no single
% condition for the merged trial.  It also contains info(1).snapshotCond(i) set
% to the condition of the trial from which the i-th snapshot was taken.  The
% returned Meta contains meta.startsOfFormerTrials, an array listing the
% snapshot indices in the returned Data that correspond to the beginning
% snapshot of each of the original trials.
%
% Example: 
% [info2,data2,meta2] = transformIDM_mergeTrials(info,data,meta) 
%
% History
% - 8/21/02 TMM Created file.
% - 9/29/04 TMM Added rmeta.startsOfFormerTrials to returned rmeta

function [rinfo,rdata,rmeta] = transformIDM_mergeTrials(info,data,meta)

  rmeta=meta;
  rj=1;
  for t=1:1:meta.ntrials
    rmeta.startsOfFormerTrials(t)=rj;
    for j=1:1:(info(t).len)
      rdata{1}(rj,:) = data{t}(j,:);
      rinfo(1).snapshotCond(rj)=info(t).cond;
      rj=rj+1;
    end;
  end;
  
  rmeta.ntrials=1;
  rinfo(1).mint=info(1).mint;
  rinfo(1).maxt=info(end).maxt;
  rinfo(1).len=size(rdata{1},1);
  rinfo(1).cond=-1;  % cond=-1 signifies no trial condition

  if (rinfo(1).len ~= (rinfo(1).maxt - rinfo(1).mint + 1)) % sanity check
    fprintf('transformIDM_mergeTrials: WARNING info.len-1 NEQ info.maxt - info.mint!');
  end
  
  if (rinfo(1).len ~= rmeta.nsnapshots) % another sanity check
    fprintf('transformIDM_mergeTrials: WARNING info.len NEQ meta.nsnapshots!');
  end
