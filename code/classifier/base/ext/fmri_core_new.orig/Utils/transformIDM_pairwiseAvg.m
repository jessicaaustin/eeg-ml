% transformIDM_pairwiseAvg(info,data,meta)
% Returns a copy of info,data, and meta, replacing each snapshot(t) by the
% average of shapshot(t) and snapshot(t+1).  Note this shortens the length of
% each trial by one snapshot. Info and data are updated accordingly.
% Returns a cell array containing the revised info,data,meta
% Example:  
%  [info,data,meta] = transformIDM_pairwiseAvg(info,data,meta);
%  
%
% Dependencies: none
%
% History
% - 8/11/02 TMM Created file.
% - 9/1/02 TMM bugfix: made rdata a proper column cell array

function [rinfo,rdata,rmeta] = transformIDM_pairwiseAvg(info,data,meta)
  rdata = cell(meta.ntrials,1);  
  rinfo=info;
  %  rdata=data;  % would this make it more efficient?
  rmeta=meta;
  for t=1:1:meta.ntrials
    for j=1:1:(info(t).len - 1)
      rdata{t}(j,:) = (data{t}(j,:) + data{t}(j+1,:)) / 2.0;
    end;
    rinfo(t).len = info(t).len - 1;
  end;
  
  % should really count the new number of snapshots to be safe...
  rmeta.nsnapshots= meta.nsnapshots-meta.ntrials; 
