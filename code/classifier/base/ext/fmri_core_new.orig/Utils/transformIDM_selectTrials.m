% transformIDM_selectTrials(info,data,meta,trials) 
%
% Returns a copy of info,data,meta containing only the specified trials.  The
% input parameter 'trials' is an array listing the indices of trials to be
% included.
%
% Example:  select just trials 3 and 5
% [info2,data2,meta2] = transformIDM_selectTrials(info,data,meta,[3,5]) 
%
% History
% - 9/1/02 TMM Created file.
% - 5/9/05 indra Update mint and maxt after the selection of trials

function [rinfo,rdata,rmeta] = transformIDM_selectTrials(info,data,meta,trials)
  ntrials=length(trials);
  
  rdata = cell(ntrials,1);  
  rmeta=meta;
  for j=1:1:ntrials
    t=trials(j);  % get trial number
    rdata{j} = data{t};
    rinfo(j) = info(t);
  end;

  rmeta.ntrials=length(trials);
  rmeta.nsnapshots= sum([rinfo.len]);

  % update mint and maxt
  tStart = 1;
  for j=1:rmeta.ntrials
    tEnd = tStart + rinfo(j).len - 1;
    rinfo(j).mint = tStart;
    rinfo(j).maxt = tEnd;
    
    tStart = tEnd + 1;
  end;