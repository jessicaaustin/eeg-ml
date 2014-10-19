% transformIDM_avgTrials(info,data,meta,trials)
% 
% Returns a copy of info,data, and meta, containing a single trial which is
% the average of all trials specified by the input.  The input 'trials' is an
% array listing the indices of trials to be included in the average.  If
% trials are of different length, prints a warning message and uses the
% length of the shortest trial.
%
% Example:  
%  [info,data,meta] = transformIDM_avgTrials(info,data,meta,[3,6]);
%   returns an IDM containing one trial which is the average of trials 3 and
%   6
%  
%
% Dependencies: none
%
% History
% - 8/29/02 TMM Created file.

function [rinfo,rdata,rmeta] = transformIDM_avgTrials(info,data,meta,trials)
			     
  % calculate minimum length of the trials that will be averaged
  len=min([info(trials).len]);
  fprintf('trial lengths:');
  fprintf(' %d',[info(trials).len]);
  fprintf('\n will return avg trial of length %d', len);

  firstTrial=trials(1);
  rinfo= [info(firstTrial)];
  rinfo(1).len=len;
  
  % rdata{1} will hold sum of trials, later divide by number of trials
  rdata{1} = data{firstTrial}(1:len,:);

  for j=2:1:length(trials)
    t=trials(j);  % get trial number, then add trial to sum
    rdata{1} = rdata{1} + data{t}(1:len,:);
  end;

  rdata{1}=rdata{1} ./ length(trials);

  rmeta=meta;
  rmeta.ntrials=1;
  rmeta.nsnapshots=rinfo(1).len;

  
% brainLex phaseI trials:
% trials=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17];  
%
% brainLex phaseIIItrials
% trials=[21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38]
