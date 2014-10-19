% [info,data,meta] = transformIDM_meanTrialCond( info,data,meta )
%
% Averages all the trials of a given condition into a single one,
% handling variable length by picking the length of the shortest
% trial.
%
% In:
% - info+data+meta
%
% Out:
% - average info+data+meta, where the data structure contains
% the average trial for each condition
% - extra stdData structure, with the same structure as the average Data.
% Each point in a voxel's time series is the standard deviation of
% the corresponding point in the voxel's time series in average Data
%
% History:
% - 26 Oct 2004 - fpereira - created
%
% Notes:
% 
% Any warnings, caveats or anything else

function [avgInfo,avgData,avgMeta,stdData] = transformIDM_meanTrialCond( varargin )

DEBUG = 0;

  % varargin lets you give the function a variable number of
  % arguments, and act according to what is given. The following
  % code is an example of how to process one mandatory and one
  % optional argument (respectively arg1 and arg2)
  
  l = length(varargin);
  if l < 3
    help transformIDM_meanTrialCond; return;
  else
    info = varargin{1};
    data = varargin{2};
    meta = varargin{3};
  end

  % figure out things
  trialSequence = [info(:).cond];
  trialLens     = [info(:).len];
  nTrials       = length(trialSequence);
  posTrials     = find(trialSequence > 0);
  conditions    = unique(trialSequence(posTrials));
  nConds        = length(conditions);
  nVoxels       = size(data{1},2);
  
  % also figure out which trials belong to each condition
  for c = 1:nConds
    cond = conditions(c);
    trialsPerCond{c}    = find(trialSequence == cond);
    trialsPerCondLen{c} = trialLens(trialsPerCond{c});
    minLenPerCond(c)    = min(trialsPerCondLen{c});
    if DEBUG
      fprintf('%d\n',cond);    
      fprintf('\t');fprintf('%d ',trialsPerCond{c});fprintf('\n');
      fprintf('\t');fprintf('%d ',trialsPerCondLen{c});fprintf('\n');
    end
      
    % all greater than 0, so inverse mapping
    condToIdx(cond) = c; 
  end
  
  %
  % Average the trials
  %
  
  avgData = cell(nConds,1);
  stdData = cell(nConds,1);
  
  for c = 1:nConds
    
    avgInfo(c).cond = conditions(c);
    avgInfo(c).len  = minLenPerCond(c);
    
    % place all trials in a 3D matrix
    tmp = zeros( minLenPerCond(c), nVoxels, length(trialsPerCond{c}) );
    
    for nt = 1:length(trialsPerCond{c})
      trial       = trialsPerCond{c}(nt);      
      tmp(:,:,nt) = data{trial}( 1:minLenPerCond(c),:);
    end
   
    % compute mean and stdev trials
    avgData{c} = mean(tmp,3);
    stdData{c} = std(tmp,0,3);
  
    clear tmp;
  end
  
  avgMeta = meta;
  avgMeta.ntrials = nConds;
    
% If you want, add a sequence of tests which return 0 if all went
% well or anything else if not. Eventually we might generate test
% suites for the code automatically from these functions...

function [errorStatus] = testThis()

template = (sin(-2*pi:.01:2*pi))';

for nt = 1:10;
  info(nt).cond = 2;
  info(nt).len  = size(template,1);
  data{nt} = template + randn(size(template));;
end
meta.ntrials = 10;

[avgInfo,avgData,avgMeta,stdData] = transformIDM_meanTrialCond(info,data,meta);

