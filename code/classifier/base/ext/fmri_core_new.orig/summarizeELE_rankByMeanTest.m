% Rank the features by how likely it is that their mean value
% across features differs from 0 (as judged by a t-test).
% This code is called by _rankByNonZeroMean and _rankByPosMean.
%
% The procedure is as follows:
%
% - Find p-values for a t-test of every feature against 0,
% separately for each condition, and rank by  mean != 0 or mean > 0
%
% - Use the ranks to select features
% Conceptually, this works by taking the per condition rankings
% Cond  2   ... last condition being considered
%       f(1)    f(1)
%       f(2)    f(2)
% and  -------------- sliding a window down until the number left
% above is the # to keep.
%
% Note that features near the top of several of the rankings will
% only get selected once.
%
% In:
% - examples
% - labels
% - test type - either 'greaterThanZero' or 'differentFromZero'
%
% Out:
% - a list of features, ordered by the rank positions described above
% - the rank positions of those features
%
% Dependencies:
%
% History
% - 10 Nov 2004 - fpereira - created
%

function [sortedFeatures,sortedRankings] = summarizeELE_rankByMeanTest(varargin)
  
  l = length(varargin);
  
  if l < 2; help summarizeELE_rankByMeanTest; return;
  else
    examples = varargin{1};
    labels   = varargin{2};
    testType = 'greaterThanZero';

    if l > 2; testType = varargin{3}; end
    switch testType
     case {'greaterThanZero'}
      tail = 1;
     case {'differentFromZero'}
      tail = 0;
     otherwise
      fprintf('summarizeELE_rankByMeanTest: testType can be "greaterThanZero" or "differentFromZero"\n');return;
    end

    removeOutliers = 0;
    if l > 3; removeOutliers = varargin{4}; end
  end

  [nExamples,nFeatures] = size(examples);
  sortedLabels = unique(labels);
  nLabels      = length(sortedLabels);
  labelPvalues = zeros(nLabels,nFeatures);
  labelRanking = zeros(nLabels,nFeatures);
  
  %
  % Find p-values for a t-test of every feature against 0, separately
  % for each condition
  %
  
  for l=1:nLabels
    label     = sortedLabels(l);
    indices   = find(labels == label);
    nindices  = length(indices);
    lexamples = examples(indices,:);

    if removeOutliers
      % replaces outliers by mean of the sample w/o them
      %means = mean(examples(indices,:),1);
      %stdvs = std(examples(indices,:),0,1);
      above = prctile(lexamples,99.5);
      below = prctile(lexamples,0.5);
      aboveMask = lexamples >= repmat(above,size(lexamples));
      belowMask = lexamples <= repmat(above,size(lexamples));
      tossMask = aboveMask + belowMask;
      keepMask = ~tossMask;
    
      lexamples      = lexamples .* keepMask;
      nleft          = sum(keepMask,1);
      meanNoOutliers = sum(lexamples,1) ./ nleft;
      tossMask       = tossMask .* (ones(nindices,1)*meanNoOutliers);
      lexamples      = lexamples + tossMask;      
    end
 
    % compute p-values and rank features
    % 0.05 argument doesn't matter here
    [h,labelPvalues(l,:)]     = computeTtest(lexamples,0,0.05,tail);  
    [dummy,labelRanking(l,:)] = sort(labelPvalues(l,:));
  end
  
  %
  % Now rank
  %

  % Conceptually, this works by taking the per condition rankings
  % Cond  2   ... last condition being considered
  %       f(1)    f(1)
  %       f(2)    f(2)
  % and  -------------- sliding a window down until the number left
  % above is the # to keep.
  %
  % Note that features near the top of several of the rankings will
  % only get selected once.

  % sort rankings into the order of features
  tmp = repmat(1:nFeatures,[nLabels,1]);
  
  % labelRanking(:) - vector where the first #labels entries are
  % the features that ended up in the first place of the per label
  % rankings, the second #labels are second places, etc.

  % resort that into feature order and ...
  tmp = sortrows([labelRanking(:),tmp(:)],1);

  % ... by reshaping, get the ranking of each feature (column) over different labels (rows)
  crossLabelRankings = reshape(tmp(:,2),[nLabels,nFeatures]);
  featureBestRank    = min(crossLabelRankings,[],1);
  
  % finally, sort features by the best rank they obtained across
  % the different labels
  [sortedRankings,sortedFeatures] = sort(featureBestRank);


  
  
  
function [] = testThis()

  % mock example set
  % - first 5 features  (nonactive->active) in condition 1
  % - second 5 features (nonactive->active) in condition 2
  part1 = randn(20,5);
  part1 = part1 + ones(20,1) * [0 0.25 0.5 0.75 1];
  part2 = randn(20,5);
  part2 = part2 + ones(20,1) * [0 0.25 0.5 0.75 1];
  examples = randn(40,10);
  examples(1:20,1:5)  = part1;
  examples(21:40,6:10) = part2;
  labels   = zeros(40,1);
  labels(1:20) = 1; labels(21:40) = 2;
  
  expinfo.featureToIDM = [(1:40)',ones(40,1)];
