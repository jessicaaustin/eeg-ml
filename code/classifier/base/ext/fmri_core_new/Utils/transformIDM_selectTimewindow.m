% transformIDM_selectTimewindow(info,data,meta,snapshots) 
%
%  Returns a copy of info,data,meta containing only the specified snapshots in
% time within each trial.  The input parameter 'snapshots' is an array
% listing the indices of snapshots to be included, assuming the index of the
% first snapshot of each trial is 1.
%
% Example:  select just snapshots 3,4,5, and 7 from each trial
% [info2,data2,meta2] = transformIDM_selectTimewindow(info,data,meta,[3,4,5,7]) 
%
% History
% - 11/2/02 TMM Created file.

function [rinfo,rdata,rmeta] = transformIDM_selectTimewindow(info,data,meta,snapshots)
  ntrials= length(data);
  rdata = cell(ntrials,1);  
  rmeta=meta;
  
  for j=1:1:ntrials
    rdata{j} = data{j}(snapshots,:);
    rinfo(j) = info(j);
    rinfo(j).len = length(snapshots);
  end;
  rmeta.nsnapshots= sum([rinfo.len]);

  
