%
% Rank features by one of several methods used by CCBI
%
% In:
% - examples
% - labels
% - CCBI method to use
%        'replicability'
%                For every pair of presentations, computes the correlation of each voxel
%                over all words in the first and second presentations in the pair. Final
%                score is the average correlation over pairs.
%
%        'slope'
%                Average all the presentations images of each word. For each voxel, sort
%                its values over the 14 average word images and fit a line to them. The
%                score is the slope of the line.
%
%        'distance'
%                Similar to slope, but computes the distance of each of the 14 values
%                to their mean, and the score is the maximum distance.
%
% Out:
% - sortedFeatures - the features ranked from best to worst
% - sortedValues   - their scores, in the same order
% - featureValues  - the score at each feature, in original feature order
%
% Notes:
%
%
% Examples:
% [sortedFeatures,sortedValues,featureScores] = summarizeELE_rankByCCBImethod(examples,labels,'replicability');
%
% History:
% - 2005 Aug 17 - fpereira
%

function [sortedFeatures,sortedValues,featureScores] = summarizeELE_rankByCCBImethod(varargin)

this = 'summarizeELE_rankByCCBImethod';

%
% Process arguments
% 

DEBUG = 0;

l = length(varargin);
if l < 3; eval(sprintf('help %s;',this)); return; end

examples  = varargin{1};
labels    = varargin{2};
method    = varargin{3};
[totalLength,nFeatures] = size(examples);

%
% Figure out things and check arguments
%

sortedConds = unique(labels); nConds = length(sortedConds);
exampleBlockSize = 1;

%
% Group examples into blocks of consecutive examples with the same label.
% (intermediate representation that is useful)

for c = 1:nConds
  cond  = sortedConds(c);  
  examplesCond{c}  = find(labels == cond);
  nExamplesCond(c) = length(examplesCond{c});
  
  if rem(nExamplesCond(c),exampleBlockSize)
    fprintf('%s: error: #examples=%d in c=%d must be divisible by %d\n',nExamplesCond(c),cond,exampleBlockSize); pause; return
  else
    nBlocks(c) = nExamplesCond(c) / exampleBlockSize;
  end

  bidx = 1;
  exampleBlocks{c} = cell(nBlocks(c),1);
  for b = 1:nBlocks(c)
    exampleBlocks{c}{b} = examplesCond{c}(bidx:(bidx+exampleBlockSize-1));
    bidx = bidx + exampleBlockSize;
  end
end

if DEBUG
  for c = 1:nConds
    cond = sortedConds(c);
    fprintf('cond=%d\n\t',cond);
    fprintf('%d ',examplesCond{c});fprintf('\n');
  end
  pause    
  for c = 1:nConds
    cond = sortedConds(c);
    fprintf('cond=%d #blocks=%d\n',cond,nBlocks(c));
  end
  pause
end

if sum(diff(nBlocks))
  fprintf('%s: error: all task conditions need to have the same # of trials\n',this);
  disp(sortedConds); disp(nBlocks); pause; return;
end

nEpochs = nBlocks(1); % how many times does each condition repeat?

%
% Reorganize data by epoch
% epochData{e} has
%   example_1_cond1
%   ...
%   example_<blockSize>_cond1
%
%   example_1_cond2
%   ...
%   example_<blockSize>_cond2
%   ...
%   example_<blockSize>_cond6
%

epochLen  = nConds * exampleBlockSize;
epochData = cell(nEpochs,1);
epochMean = cell(nEpochs,1);
epochStdv = cell(nEpochs,1); 
avgEpoch  = zeros(epochLen,nFeatures);

% Organize trials per epoch
% (each epoch contains 1 trial of each condition)
% and discard original data to save memory

for e = 1:nEpochs

  epochData{e} = zeros(epochLen,nFeatures);

  % orders examples by condition
  idx = 1;
  for c = 1:nConds
    indices = exampleBlocks{c}{e};
    epochData{e}(idx:(idx+exampleBlockSize-1),:) = examples(indices,:); 
    idx = idx + exampleBlockSize;
  end

  epochMean{e} = mean(epochData{e},1);
  epochStdv{e} = std( epochData{e},0,1);
  avgEpoch     = avgEpoch + epochData{e};
end
clear examples;

avgEpoch = avgEpoch / nEpochs;

%
% Compute rankings
%

switch method
  
  
 case {'slope','tbetaslope'}

  % same, except one uses a single average epoch and the other
  % all the epochs
  
  if isequal(method,'slope')
    tmpData{1} = avgEpoch;
  else
    tmpData    = epochData;
  end
  nToUse = length(tmpData);

  epochIntercepts = zeros(nToUse,nFeatures);
  epochSlopes     = zeros(nToUse,nFeatures);
  
  for t = 1:length(tmpData)
  
    tmp = sort(tmpData{t},1);
    len = size(tmpData{t},1);
  
    x = [ones(len,1),(1:len)'];

    fprintf('%s: computing slopes\t',this);
    for f = 1:nFeatures  
      %tmp
      %pause
      %      if rem(f,100); fprintf('%d ',f); end
      coefficients = regress(tmp(:,f),x);
      epochIntercepts(t,f) = coefficients(1);
      epochSlopes(t,f)     = coefficients(2);
    end
    fprintf('\n');
  end

  if isequal(method,'slope')
    featureScores = epochSlopes; % only one line
  else
    % t-ttest on whether the mean of epochSlopes is > 0
    % NOT IMPLEMENTED YET, use mean only
    featureScores = mean(epochSlopes,1);
  end

  
 case {'distance'}
  
  tmp = mean(avgEpoch,1);
  len = size(avgEpoch,1);
  featureMaxDistFromMean = max( abs( avgEpoch - repmat(tmp,len,1)), [],1 );
  featureScores = featureMaxDistFromMean;
 
  
 case {'replicability'}
%keyboard
  similarityMeasure = 'correlation';
  
  % make each epoch mean image 0  
  for e = 1:nEpochs; epochData{e} = epochData{e} - repmat(epochMean{e},epochLen,1);end

  %% Compute mean correlation of each voxel with itself in each pair of epochs
  
  nPairs = (nEpochs * (nEpochs-1)) / 2;
  epochPairVoxelCorrelations = zeros(nEpochs,nFeatures);

  pidx = 1;
  for e1 = 1:(nEpochs-1)
    for e2 = (e1+1):nEpochs
      
      switch similarityMeasure
       case {'correlation','correlationWithMean'}
	% compute correlation for each voxel between the two epochs in this pair    
	epochPairVoxelCorrelations(pidx,:) = sum(epochData{e1} .* epochData{e2}) ./ ((epochStdv{e1}.*epochStdv{e2})*(epochLen-1));
       case {'euclidean'}
	epochPairVoxelCorrelations(pidx,:) = sqrt(sum((epochData{e1}-epochData{e2}).^2,1));
       case {'absolute'}
	epochPairVoxelCorrelations(pidx,:) = sqrt(sum(abs(epochData{e1}-epochData{e2}),1));
      end
      
      pidx = pidx + 1;
    end
  end

  meanPairVoxelCorrelations = mean(epochPairVoxelCorrelations,1);
  featureScores = meanPairVoxelCorrelations;
  
 otherwise
  
end

%
% Postprocess rankings so that feature rankings are from best to worst
%
  

switch method

 case {'slope','tbetaslope','distance','replicability'}
  % large is good
  [sortedValues,sortedFeatures] = sort(featureScores); % small first
  sortedFeatures = fliplr(sortedFeatures); % large first
  sortedValues   = fliplr(sortedValues);
  
 otherwise
  [sortedValues,sortedFeatures] = sort(featureScores); % small first
  
end
  

return;






function [] = testThis()

% example from the IDM test, below is code
% to transform the IDM into examples

stamp = [ [1 2 3 4 5 6 7 8]',...
	  [8 7 6 5 4 3 2 1]',...
	  [4 4 4 4 4 4 4 4]'];

zstamp = zeros(size(stamp));



nFeatures = 12;
nEpochs = 3;
len     = 8;

nt = 1;

% first six voxels are consistent, next six aren't

info(nt).cond = 2;
info(nt).len  = len;
data{nt} = [stamp,zstamp,stamp,zstamp];
data{nt} = data{nt} + randn(size(data{nt}));
nt = nt + 1;

info(nt).cond = 3;
info(nt).len  = len;
data{nt} = [zstamp,stamp,zstamp,stamp];
data{nt} = data{nt} + randn(size(data{nt}));
nt = nt + 1;

info(nt).cond = 3;
info(nt).len  = len;
data{nt} = [zstamp,stamp,zstamp,zstamp];
data{nt} = data{nt} + randn(size(data{nt}));
nt = nt + 1;

info(nt).cond = 2;
info(nt).len  = len;
data{nt} = [stamp,zstamp,zstamp,zstamp];
data{nt} = data{nt} + randn(size(data{nt}));
nt = nt + 1;

info(nt).cond = 2;
info(nt).len  = len;
data{nt} = [stamp,zstamp,zstamp,zstamp];
data{nt} = data{nt} + randn(size(data{nt}));
nt = nt + 1;

info(nt).cond = 3;
info(nt).len  = len;
data{nt} = [zstamp,stamp,zstamp,zstamp];
data{nt} = data{nt} + randn(size(data{nt}));
nt = nt + 1;


% transform this into examples by concatenating all the trials
meta.ntrials = 6;
[dummy,d,dummy] = transformIDM_joinTrials(info,data,meta);
examples = d{1};
labels = [];
for nt = 1:meta.ntrials
  labels = [labels;repmat(info(nt).cond,info(nt).len,1)];  
end

[sortedFeatures,sortedValues] = summarizeELE_rankByConsistency(examples,labels)

% now use this to check voxels 1-6 by hand

v = 3; p1 = [data{1}(:,v);data{2}(:,v)]; p2 = [data{5}(:,v);data{6}(:,v)];corrcoef(p1,p2)

