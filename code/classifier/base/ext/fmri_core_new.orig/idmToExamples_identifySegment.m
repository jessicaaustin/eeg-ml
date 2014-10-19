function [examples, labels, expInfo] = idmToExamples_identifySegment(varargin)
% Given an IDM, create examples to be identified as a picture or 
% sentence segment
%
% In:
% - IDM and an optional argument sepecifying whether to normalize the data
%   The default setting is 'normalize'
%
% Out:
% - examples, labels and experimental information
%
% Dependencies:
%
% Examples:
% - [examples, labels, expInfo] = idmToExamples_identifySegment(info, data, meta, 0)
%
% History:
% - 18 Oct. 2002 Xuerui Created from exp_starplus_identifySegment.m in
%                       fmri_code/Classification

l = length(varargin);
normalize = 1;
if l < 3
    fprintf(1,'syntax: idmToExamples_identifySegment(info,data,meta,<normalize>)');
    examples = []; labels = [];
    return;
else
    info = varargin{1};
    data = varargin{2};
    meta = varargin{3};
    study = meta.study;
end

if l > 3
    normalize = varargin{4};
end

expInfo.experiment = 'identifySegment';
expInfo.normalized = normalize;
expInfo.meta = meta;

% gather information about number,type and length of trials
%[ntrials,nvoxels,nconds,minTrialLenCond,ntrialsCond] = mri_infoTrials(info, data, meta, study);
IDM_information=IDMinformation( info,data,meta,meta.study );
ntrials=IDM_information.nTrials;
nvoxels=IDM_information.nVoxels;
nconds=IDM_information.nConds;
minTrialLenCond=IDM_information.minTrialLenCond;
ntrialsCond=IDM_information.nTrialsCond;
  
%if normalize
%    [info, data, meta] = transformIDM_normalizeInterval(info, data, meta, 0, [1 32]);
%end
 if normalize
     % normalize to [0,1]
     for i = 1:1:ntrials
         if info(i).cond > 1
             for j = 1:1: nvoxels
                 mn = min(data{i}(1:32,j));
                 mx = max(data{i}(1:32,j));
                 data{i}(1:32,j) = (data{i}(1:32,j)-mn)./(mx-mn);
             end
         end
     end  
 end

minTrialLen = min(minTrialLenCond(2:1:nconds));

%
% set up for outputting examples
%
% Examples are stored as 
% v1...vn (time1) ... v1 ... vn (time t)
%

nfeatures = 16 * nvoxels; 
nexamples = 2 * (ntrialsCond(2) + ntrialsCond(3));
nlabels   = 1;

examples = zeros(nexamples,nfeatures);
labels   = zeros(nexamples,nlabels);

% go ahead and do it
exampleIdx = 1;
for nt = 1:1:ntrials
    len  = info(nt).len;
    cond = info(nt).cond;
    
    if cond > 1
        
        % Create the example for the first 16 time points
        featureIdx = 1;
        % concatenate the signal at all voxels for each time point
        for t = 1:1:16
            featureIdxNext = featureIdx + nvoxels;
            examples( exampleIdx, featureIdx:1:(featureIdxNext-1) ) = data{nt}(t,:);
            featureIdx = featureIdxNext;
        end
        if isequal(study,'data-starplus-sp')
            labels(exampleIdx,1) = 1;
        else
            labels(exampleIdx,1) = 0;
        end
        exampleIdx = exampleIdx + 1;
        
        % Create the example for the next 16 time points
        featureIdx = 1;
        % concatenate the signal at all voxels for each time point
        for t = 17:1:32
            featureIdxNext = featureIdx + nvoxels;
            examples( exampleIdx, featureIdx:1:(featureIdxNext-1) ) = data{nt}(t,:);
            featureIdx = featureIdxNext;
        end      
        if isequal(study,'data-starplus-sp')
            labels(exampleIdx,1) = 0;
        else
            labels(exampleIdx,1) = 1;
        end
        exampleIdx = exampleIdx + 1;        
    end
end













