% transformIDM_unfold(info,data,meta)
% Returns an IDM where each trial gets turned into |trial| one-time point
% trials with the same condition. E.g. a trial with condition 3 and length
% 20 becomes a succession of 20 trials with condition 3 and length 1.
% Info,data and meta are updated accordingly.
% This is used in block studies to turn each time point in a block into
% a trial, and then each trial into an example.
% Example:  
%  [info,data,meta] = transformIDM_unfold(info,data,meta); 
%
% Dependencies:
% - transformIDM_unfold depends on mri_infoTrials
%
% History
% - 10 Aug 02 - fpereira - created from process_dataToExamples
%
%

function [uinfo,udata,umeta] = transformIDM_unfold(info,data,meta)
  
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
   
  
  udata=[];
  umeta=[];
  
  idxTrial = 1;
  for nt=1:ntrials
      len=info(nt).len;
      cond=info(nt).cond;
      for inTrialIdx = 1:len
         currTrialInfo=info(nt);
         currTrialInfo.len=1;
         if isfield(currTrialInfo,'mint')
             currTrialInfo.mint = idxTrial;
         end
         if isfield(currTrialInfo,'maxt')
             currTrialInfo.maxt = idxTrial;
         end
         if isfield(currTrialInfo,'nmissing')
             currTrialInfo.nmissing = 0;
         end
         uinfo(idxTrial)=currTrialInfo;
         
         udata{idxTrial}=data{nt}(inTrialIdx,:);
         
         idxTrial=idxTrial+1;
      end
  end
  
  umeta=meta;
  umeta.ntrials=idxTrial-1;
  
  % info+data+meta were transformed into a pseudo-dataset
  % with as many trials as images (minus the ones excluded)
  % and where each trial has length 1   
  
