% K-nearest neighbours
%
% Allows specification of several different distance metrics and
% also variable weighing of dimensions in a few of those.
%
% In:
% - k
% - training examples + labels
% - test examples
% - (optionally) distance metric - defaults to euclidean distance
%
% Out:
% - A score matrix with dimensions #test examples x #labels.
%   cell(e,j) has the fraction of the knneighbours of e that
%   have the j-th label (in sorted order)
%
% Dependencies:
%
%
% History
% - 11 Mar 05 - fpereira - adapted from a previous version
%
%

function [Scores] = classifierKNN2( varargin )

%
% Process parameters
%

l = length(varargin);
if l < 4; help classifierKNN2; return;
else
  k = varargin{1};
  TrainExamples = varargin{2};
  TrainLabels   = varargin{3};
  TestExamples  = varargin{4};
  [numTrain,numFeatures] = size(TrainExamples);
  numTest = size(TestExamples,1);
  
  if l > 4; additionalParameters = 1; else; additionalParameters = 0; end
end

% figure out a few things
sortedLabelValues = sort(unique(TrainLabels));
numLabels         = length(sortedLabelValues);
  
% Compute distances between all test examples and training examples
if ~additionalParameters
  [Distances] = compute_exampleDist(TestExamples,TrainExamples);
else
  % build a command line that passes the additional parameters
  cmd = 'TestExamples,TrainExamples';
  for m=5:l; cmd = sprintf('%s,varargin{%d}',cmd,m); end
  cmd = sprintf('[Distances] = compute_exampleDist(%s);',cmd);
  eval(cmd);
end

% Need to know how many of each label are close to a test example

% For each test example, sort indices of train examples by distance to it
[SDistances,SIndices] = sort(Distances,2);
Scores = zeros(numTest,numLabels);
  
% Now create scores over target labels. Normally this
% is just the fraction of k-nearest neighbours that
% have a given label

if ~isequal(sort(TrainLabels),sortedLabelValues)
  
  IndicesKNN = SIndices(:,1:k); % indices of knn
  % And find their labels
  tmp = TrainLabels(IndicesKNN(:)); % all labels
  LabelsKNN  = reshape(tmp,size(IndicesKNN));
  
  for l=1:numLabels
    label = sortedLabelValues(l);
    % counts the number of times each label appears in the knn of
    % each example
    tmp = (LabelsKNN == label);
    Scores(:,l) = sum(tmp,2);
  end
  
  Scores = Scores ./ k; % make the number a fraction
else
  % but when having a single example from every class
  % in the training set, the actual ranking by distance
  % is important, so here the scores are actually -1*distances
     
  tmp = [Distances',TrainLabels];
  tmp = sortrows(tmp,[numTest+1]); % order by label
  Scores = -1 * (tmp(:,1:numTest)');
end



  
  
function [] = testThis()
  
  trainSet    = zeros(5,2);
  trainLabels = zeros(5,1);
  
  trainSet(1,:) = [1 2];
  trainSet(2,:) = [1 4];
  trainSet(3,:) = [3 1];
  trainSet(4,:) = [3 3];
  trainSet(5,:) = [3 5];
  trainLabels(:,1) = ([0 0 1 1 1])';
  ntrain = size(trainSet,1);
  class0idx = find(trainLabels == 0);
  class1idx = find(trainLabels == 1);
  plot(trainSet(:,1),trainSet(:,2),'.','MarkerSize',1);
  axis([0 6 0 6]);
  text(trainSet(class0idx,1),trainSet(class0idx,2),'-');
  text(trainSet(class1idx,1),trainSet(class1idx,2),'+');  
  
  
  testSet = zeros(2,2);
  testSet(1,:) = [2 3];
  testSet(2,:) = [5 3];
  testLabels  = ([0 1])';
  ntest  = size(testSet,1);

  for t=1:1:ntest
    fprintf(1,'testing: %d\t%d\t',testSet(t,1),testSet(t,2));
    
    for k=1:1:1
      [testLabels] = classifier_knn( k, trainSet, trainLabels, testSet);
    end
    
    fprintf(1,'\n');
  end
