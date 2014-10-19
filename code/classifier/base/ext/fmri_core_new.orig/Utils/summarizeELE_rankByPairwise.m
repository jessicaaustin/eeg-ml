% Take each pair of classes and rank the features by how well they
% do at discriminating the examples belonging to that pair.
% Then combine those rankings into a single ranking in the
% following way. Take the top feature in each pairwise ranking,
% rank these first. Then take the second best from each ranking,
% third best, etc.
%
% This uses summarizeELE_rankByNestedCV to rank how well the features
% distinguish each pair of classes, please read the help for that.
%
% In:
% - examples (#examples x #features array)
% - labels 
% - optional
%   - classifier to use (only nbayes is supported right now,defaults to that)
%   - error measure to use - defaults to 'error', could also be 'averageRank'
%
% Out:
% - sortedFeatures - features sorted from best to worst
% - sortedRankings   - feature ranking position sorted from best to worst (matches sortedFeatures)
%
% Examples:
% - [sortedFeatures,sortedRankings] = summarizeELE_rankByPairwise( examples,labels,'nbayes' )
%
% History
% - 2005 Aug 31 - fpereira - created from older code

function [sortedFeatures,sortedRankings] = summarizeELE_rankByPairwise(varargin)
  
  l = length(varargin);
  
  if l < 2; help summarizeELE_rankByPairwise; return;
  else
    examples   = varargin{1};
    labels     = varargin{2};
    classifier = 'nbayes';
    classifierParameters = {};
    errorMeasure = 'error';
    
    if l > 3
      classifier = varargin{3};
      if l > 4
	errorMeasure = varargin{4};
      end
    end
  end
      
  % figure out a few things
  sortedLabelValues = sort(unique(labels)); nLabels = length(sortedLabelValues);
  nPairs = (nLabels*(nLabels-1))/2;
  
  [nExamples,nFeatures] = size(examples);
  
  pairRanking = zeros(nPairs,nFeatures);
  pairScores  = zeros(nPairs,nFeatures);
  
  %
  % Compute nested cross validation error for every feature on the
  % pairwise discrimination task
  %
  
  idx = 1;
  for l1 = 1:nLabels-1
    
    label1 = sortedLabelValues(l1);
    
    for l2 = l1+1:nLabels
      
      label2 = sortedLabelValues(l2);
      pairLabelValues = [label1 label2];      
      fprintf('processing %d and %d\n',label1,label2);

      %% find examples with these labels
      indicesTrain  = find((labels  == label1) | (labels == label2));
      
      pairExamplesTrain  = examples(indicesTrain,:);
      pairLabelsTrain    = labels(indicesTrain);
      
      %% rank them
      [sortedFeatures,sortedMeasure] = summarizeELE_rankByNestedCV(pairExamplesTrain,pairLabelsTrain,classifier,errorMeasure);

      %% resort into feature order
      pairRanking(idx,:)  = sortedFeatures;
      [dummy,resortOrder] = sort(sortedFeatures);
      pairScores(idx,:)   = sortedMeasure(resortOrder);
      
      idx = idx + 1;
    end; % over label 2    
  end; % over label 1
    
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

  % Each line of pairRankings contains the ranking of features for
  % tests involving that pair
  
  % sort those rankings into the order of features
  % tmp is:
  % 1       1 2 3 ... # features
  % 2       1 2 3 ... # features
  % ...
  % #pairs  1 2 3 ... # features  
  tmp = repmat(1:nFeatures,[nPairs,1]);
  
  % pairRanking(:) - vector where the first #labels entries are
  % the features that ended up in the first place of the per label
  % rankings, the second #labels are second places, etc.

  % resort that into feature order and ...
  tmp = sortrows([pairRanking(:),tmp(:)],1);

  % ... by reshaping into #pairs x #features, get the ranking of
  % each feature (column) over different labels (rows)
  crossLabelRankings = reshape(tmp(:,2),[nPairs,nFeatures]);
  featureBestRank    = min(crossLabelRankings,[],1);
  
  % finally, sort features by the best rank they obtained across
  % the different labels
  [sortedRankings,sortedFeatures] = sort(featureBestRank);


  
function [] = testRankingMechanism

% 3 labels, 3 pairs, 5 features

pairRanking = [
     1     3     5     2     4
     5     4     3     2     1
     1     3     2     4     5];

tmp = repmat(1:nFeatures,[nPairs,1]);
tmp = sortrows([pairRanking(:),tmp(:)],1);

% give us a #features x 2 array, organized by feature
% and containing, for each feature, its position on
% rankings for all the pairs
% e.g.
% tmp = [
%     1     1
%     1     1
%     1     5
%     2     3
%     2     4
%     ...
%     5     5]
% ... by reshaping into #pairs x #features, get the ranking of
% each feature (column) over different labels (rows)

crossLabelRankings = reshape(tmp(:,2),[nPairs,nFeatures]);
featureBestRank    = min(crossLabelRankings,[],1);
  
  
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
