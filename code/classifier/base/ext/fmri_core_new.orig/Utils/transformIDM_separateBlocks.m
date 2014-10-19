% transformIDM_separateBlocks(info,data,meta)
% Returns an IDM where each trial with a given condition number acquires
% a different condition number. For instance, the sequence of conditions
% in a study with 6 conditions
% 2,3,2,4,3,5,6,3
% would become
% 2,3,7,4,8,5,6,9
% Info is updated accordingly.
% One application of this was to transform 6 into 12 categories in the
% data from the Categories study.
% Example:
%   [info,data,meta] = transformIDM_separateBlocks(info,data,meta);
%
%
% WARNING:
% Note that it has been designed with the Categories study in mind,
% where each there are two trials in each condition, each of which
% is a block. If there aren't two trials (subjects 05695 and
% 05675), this will still work but the condition numbers won't
% match those in other subjects.
%
% WARNING: actually, now it assumes an even number of trials so
% that we can keep track of what each trial/condition get turned into
%
% Dependencies:
% - transformIDM_separateBlocks depends on mri_infoTrials
% - transformIDM_separateBlocks depends on
%
% History
% - 13 Aug 02 - fpereira - created from process_dataToExamples
%
%

function [ninfo,data,meta] = transformIDM_separateBlocks(info,data,meta)
  
% gather information about number,type and length of trials and other dataset parameters
  IDM_information=IDMinformation( info,data,meta,meta.study );
  ntrials=IDM_information.nTrials;
  nvoxels=IDM_information.nVoxels;
  nconds=IDM_information.nConds;
  minTrialLenCond=IDM_information.minTrialLenCond;
  ntrialsCond=IDM_information.nTrialsCond;
  
  information = localInformation(meta.study);
  dimensions  = information.dimensions;
  dimx        = dimensions(1); dimy = dimensions(2); dimz = dimensions(3);
  % check how many trials per condition and whether all conditions
  % have the same number

  a = sum(ntrialsCond(2:nconds) == ntrialsCond(2));
  if (a ~= (nconds-1))
    fprintf('transformIDM_separateBlocks: ERROR: all conditions must have the same number of trials\n');
    fprintf('\tright now:\n');
    for c=1:nconds; fprintf('\t%d\t%d\n',c,ntrialsCond(c)); end;
    return;
  end
  
  ntc       = ntrialsCond(2); % # trials per condition (not fixation)
  ncondsNew = ntc*(nconds-1) + 1; % new # conditionsQ

  % New classes will be numbered by
  % nconds + #trial(within cond)+ncond
  % 1 -> 1
  % first 2  -> 2
  % first 7  -> 7
  % second 2 -> 7 + 2 - 1  = 8
  % second 7 -> 7 + 7 - 1  = 13
  % ...
  % seventh 2 -> (6)*7 + 2 - 6 = 38 
  % seventh 7 -> (6)*7 + 7 - 6 = 43
  % eighth  2 -> (7)*7 + 2 - 7 = 44
  % eighth  2 -> (7)*7 + 7 - 7 = 49
  % ninth   2 -> (8)*7 + 2 - 8 = 50
  % so
  % twc = #trial(within,cond)
  % twc,condition => (twc-1)*#conds + condition - (twc-1)  

  ninfo   = info; % will store new condition information
  twc = zeros(nconds,1);
  
  fprintf('transformIDM_separateBlocks:\n');
  % now go through the trials and modify info accordingly
  for nt=1:ntrials
    cond = info(nt).cond;
    len  = info(nt).len;
    
    if cond < 2
      % do nothing, leave ninfo(nt) as it is
    else
      twc(cond)      = twc(cond) + 1;    
      ninfo(nt).cond = nconds*(twc(cond)-1) + info(nt).cond - (twc(cond)-1);
    end

    r = mod(ninfo(nt).cond,nconds);
    
    fprintf('\ttrial %d:\t%d\t->\t%d\n',nt,info(nt).cond,ninfo(nt).cond);
  end
  
  
