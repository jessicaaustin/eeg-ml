% Ranks features by using each one as a classifier in a
% cross-validation within the training set
% 
% In:
% - examples (#examples x #features array)
% - labels 
% - classifier to use (only nbayes is supported right now)
% - optional
%   - error measure to use - defaults to 'error', could also be 'averageRank'
%
% Out:
% - sortedFeatures - features sorted from best to worst
% - sortedErrors   - feature scores sorted from best to worst (matches sortedFeatures)
% - measurePerClass - not documented yet
%
% Examples:
% - [sortedFeatures,sortedErrors,measurePerClass] = summarizeELE_rankByNestedCV( examples,labels,'nbayes' )
%
% History
% - 2005 Aug 31 - fpereira - created from older code
%

function [sortedFeatures,sortedErrors,measurePerClass] = summarizeELE_rankByNestedCV( varargin )

%% process arguments and figure out/test a few things

l = length(varargin);
if l < 3; help summarizeELE_rankByNestedCV; return; end

examples   = varargin{1};
labels     = varargin{2};
classifier = varargin{3};
classifierParams = {};
switch classifier
 case {'nbayes','nbayesPooled'}
  % OK
 otherwise
  fprintf('summarizeELE_rankByNestedCV: classifier %s is not supported\n');return;
end
measure = 'error';
if l > 3; measure = varargin{4}; end

doFeatureSelection = 0;

%if l > 3
%  classifierParams = varargin{4};
%  if l > 4
%    doFeatureSelection = 1;
%    FSinfo   = varargin{5};
%    FSmethod = FSinfo{1};
%    FSparams = FSinfo{2};    
%  end
%end

sortedLabelValues = unique(labels);
nLabels = length(sortedLabelValues);
[nExamples,nFeatures] = size(examples);


%% Separate examples into classes and compute leave-1-out means/stdevs

indicesLabel  = cell(nLabels,1);
examplesLabel = cell(nLabels,1);
nPerLabel     = zeros(nLabels,1);

for l = 1:nLabels  
  % find examples with this label
  label = sortedLabelValues(l);  
  indicesLabel{l} = find(labels == label);
  nPerLabel(l)    = length(indicesLabel{l});
end

% ensure we have a balanced dataset
if sum(diff(nPerLabel))
  fprintf('error: there must be the same number of examples in each class\n');return;
end
nPerLabel = nPerLabel(1);

%
% Compute matrix of means without each of the points
% 

for l = 1:nLabels  
  % extract examples with this label into a separate matrix
  examplesLabel{l} = examples(indicesLabel{l},:);
  
  % get means and stdevs for each feature, leaving out each of
  % the examples (e.g. row <r> contains mean/stdev of all features,
  % over all examples except the <r>-th)
  matrixOfMeans{l}  = aux_meansMinusOneRow( examplesLabel{l} );
  matrixOfStdevs{l} = aux_stdevsMinusOneRow( examplesLabel{l} );

  % for now, hardwire class priors
  matrixOfPriors{l} = ones(nPerLabel,nFeatures);
end

%
% Compute matrix of standard deviations
%

switch classifier
  
 case {'nbayes'}
  % all set

 case {'nbayesPooled'}
  % the matrix of standard devs shouls be computed over all points
  % minus their class mean
  
end


clear examples;


%% Now compute individual feature scores

fs = cell(nLabels);

for lB = 1:nLabels

  % for the examples with label <lB>, compute log(prob(class<lA>|feature)),
  % for all classes <lA>
  
  fs{lB} = zeros(nPerLabel,nFeatures,nLabels);

  % using the MLEs for each feature and each training set in class <lA>
  % compute the feature scores    

  for lA = 1:nLabels
 
    tmp = examplesLabel{lB};
    tmp = tmp -  matrixOfMeans{lA};
    tmp = tmp ./ matrixOfStdevs{lA};
    tmp = -1 * (tmp .^2) / 2;
    tmp = tmp - log(sqrt(2*pi)*matrixOfStdevs{lA}) + log(matrixOfPriors{lA});
  
    fs{lB}(:,:,lA) = tmp;  
  end
  
end
  

%% And label the examples and compute errors

measureOverall  = zeros(1,nFeatures);
measurePerClass = cell(nLabels,1);

for lB = 1:nLabels
  
  label = sortedLabelValues(lB);
  
  % for examples of class <lB>, sort label scores
  
  % last slice of sortedIndices contains z-coordinates of highest scores
  [sortedScores,sortedIndices] = sort(fs{lB},3); % ascending
  sortedScores  = flipdim(sortedScores,3); % now descending
  sortedIndices = flipdim(sortedIndices,3);
  labelRankings   = sortedLabelValues( sortedIndices );
  predictedLabels = labelRankings(:,:,1);

  switch measure
   case {'error'}
    measurePerClass{lB} = (predictedLabels ~= label);
   case {'averageRank'}
    % labelRankings has, for each example (x) and feature (y)
    % a ranking (z) of the labels by score (highest is first)
    
    % Same dimensions, 0 everywhere except 1 where the true labels are
    % in the ranking. The position of the 1 in each ranking is the position
    % of the correct label in that ranking
    labelRankingsMask = (labelRankings == label);

    % finds the position of the 1
    [dummy,correctLabelRank] = max(labelRankingsMask,[],3);
    
    % compute normalized rank (varies between 0 (perfect) and 1)
    measurePerClass{lB} = (correctLabelRank - 1) / (nLabels - 1);
  end

  measureOverall = measureOverall + sum(measurePerClass{lB},1);
end

% average measure over all examples
measureOverall = measureOverall / nExamples;

[sortedErrors,sortedFeatures] = sort(measureOverall);


%
% Takes a matrix and computes mean across all rows except one.
% The computation is done for every row in turn (see line below
% for details)
%
% In:
% - <m> x <n> matrix to perform the computation on
%
% Out:
% - <m> x <n> matrix, where row <r> contains the mean over all rows
%   except row <r>
%
 
function [matrixOfMeans] = aux_meansMinusOneRow( m )
 
[nrows,ncols] = size(m);
matrixOfMeans = zeros(nrows,ncols);
 
matrixSum = sum(m,1);
 
for r = 1:nrows
  matrixOfMeans(r,:) = matrixSum - m(r,:);
end
 
matrixOfMeans = matrixOfMeans ./ (nrows-1);


%
% Takes a matrix and computes std across all rows except one.
% The computation is done for every row in turn (see line below
% for details)
%
% In:
% - <m> x <n> matrix to perform the computation on
%
% Out:
% - <m> x <n> matrix, where row <r> contains the standard deviation
%   over all rows, except row <r>
%
                                                                                                 
function [matrixOfStdevs] = aux_stdevsMinusOneRow( m )
                                                                                                 
[nrows,ncols] = size(m);
matrixOfStdevs = zeros(nrows,ncols);
                                                                                                 
m1 = m;
m2 = m .^ 2;
                                                                                                 
matrixOfMeans1 = aux_meansMinusOneRow(m1);
matrixOfMeans2 = aux_meansMinusOneRow(m2);
                                                                                                 
matrixOfStdevs = matrixOfMeans2 - (matrixOfMeans1 .^ 2);
matrixOfStdevs = sqrt(matrixOfStdevs);





function [] = testThis() 

%
% Simple two class example, where 1/3 of the features have small
% error, 1/3 medium and 1/3 large
% 

% generate dataset of features with 0.5 bayes error
nExamplesPerClass = 100;
nFeatures         = 300;
meanDif = 2;

examples = zeros(nExamplesPerClass*2,nFeatures);

% 1/3 small error
sidx = 1; eidx = sidx + nFeatures/3 - 1;
examples(1:nExamplesPerClass,sidx:eidx) = randn(nExamplesPerClass,nFeatures/3);
examples(nExamplesPerClass+1:end,sidx:eidx) = randn(nExamplesPerClass,nFeatures/3) + 3;

% 1/3 medium error
sidx = eidx+1; eidx = sidx + nFeatures/3 - 1;
examples(1:nExamplesPerClass,sidx:eidx) = randn(nExamplesPerClass,nFeatures/3);
examples(nExamplesPerClass+1:end,sidx:eidx) = randn(nExamplesPerClass,nFeatures/3) + 1;

% 1/3 large error
sidx = eidx+1; eidx = sidx + nFeatures/3 - 1;
examples(1:nExamplesPerClass,sidx:eidx) = randn(nExamplesPerClass,nFeatures/3);
examples(nExamplesPerClass+1:end,sidx:eidx) = randn(nExamplesPerClass,nFeatures/3) + 0.001;

labels   = ones(nExamplesPerClass*2,1);
labels(nExamplesPerClass+1:end,1) = 2;

errorMeasure = 'averageRank';
%errorMeasure = 'error';
[sortedFeatures,sortedValues] = summarizeELE_rankByNestedCV(examples,labels,'nbayes',errorMeasure);
[dummy,resorted] = sort(sortedFeatures);
plot(sortedValues(resorted));


%
% Same, but now we have three classes
% 

% generate dataset of features with 0.5 bayes error
nExamplesPerClass = 100;
nFeatures         = 300;
meanDif = 2;

examples = zeros(nExamplesPerClass*3,nFeatures);

% 1/3 small error
sidx = 1; eidx = sidx + nFeatures/3 - 1;
examples(1:nExamplesPerClass,sidx:eidx)                     = randn(nExamplesPerClass,nFeatures/3);
examples(nExamplesPerClass+1:2*nExamplesPerClass,sidx:eidx) = randn(nExamplesPerClass,nFeatures/3) + 3;
examples(2*nExamplesPerClass+1:end,sidx:eidx)               = randn(nExamplesPerClass,nFeatures/3) + 6;

% 1/3 medium error
sidx = eidx+1; eidx = sidx + nFeatures/3 - 1;
examples(1:nExamplesPerClass,sidx:eidx)                     = randn(nExamplesPerClass,nFeatures/3);
examples(nExamplesPerClass+1:2*nExamplesPerClass,sidx:eidx) = randn(nExamplesPerClass,nFeatures/3) + 1;
examples(2*nExamplesPerClass+1:end,sidx:eidx)               = randn(nExamplesPerClass,nFeatures/3) + 2;

% 1/3 large error
sidx = eidx+1; eidx = sidx + nFeatures/3 - 1;
examples(1:nExamplesPerClass,sidx:eidx)                     = randn(nExamplesPerClass,nFeatures/3);
examples(nExamplesPerClass+1:2*nExamplesPerClass,sidx:eidx) = randn(nExamplesPerClass,nFeatures/3) + 0.05;
examples(2*nExamplesPerClass+1:end,sidx:eidx)               = randn(nExamplesPerClass,nFeatures/3) + 0.1;


labels   = ones(nExamplesPerClass*3,1)*2;
labels(nExamplesPerClass+1:2*nExamplesPerClass,1) = 3;
labels(2*nExamplesPerClass+1:3*nExamplesPerClass,1) = 4;

errorMeasure = 'averageRank';
%errorMeasure = 'error';
[sortedFeatures,sortedValues] = summarizeELE_rankByNestedCV(examples,labels,'nbayes',errorMeasure);
[dummy,resorted] = sort(sortedFeatures);
plot(sortedValues(resorted));




