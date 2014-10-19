% Given scores over labels for test examples, use those scores
% to predict labels and evaluate the predictions using some
% measure, such as  error, average rank, precision/recall
% or any other.
%
% In:
% - test examples
% - models produced  by trainClassifier
% - scores structure produced by applyClassifier
% - classifier used
% - error measure 
% - test labels
%
% Out:
% - a result cell array
%
%   The first cell always contains the measure that is to be computed
%   (for accuracy and averageRank, this is a single number, but
%   could be some other structure for something like a precision
%   recall curve, for instance).
%
%   The second one contains the measure value per example, if that
%   can be calculated for the given measure (i.e. for accuracy it
%   is 1 or 0 depending on whether the classification is right or
%   wrong, for averageRank it is a number in [0,1] corresponding to
%   the rank of the correct label for that example, etc).
%
% - predictedLabels - a vector with the labels learnt
%
% - trace - a cell array containing a number of trace structures:
%
%   - labelRankings - matrix with #examples rows and #classes columns.
%     Each row contains, for one example, all the label values,
%     ranked by their decreasing scores. The most probable/highest
%     scoring label is thus that on the first column.
%
%   - labelCounts - matrix with dimensions #classes*#classes.
%     It's a confusion matrix, i.e. counts how often examples with
%     true label <row> got classified as <col>. The row and column
%     numbers can be used as indices into sortedLabelValues to
%     get the label values.
%
% Dependences:
%
% Notes:
%

% History
% - 10 Mar 2005 - fpereira - adapted from previous version
%

function [result,predictedLabels,trace] = summarizePredictions( varargin )

%
% Process parameters
%

l = length(varargin);
if l < 3; help summarizePredictions; return;
else

  scores            = varargin{1};
  trainedClassifier = varargin{2};
  measure           = varargin{3};

  if l > 3
    % right now, all the measures required true labels
    testLabels = varargin{4};
  else
    fprintf('summarizePredictions: ERROR: please provide the true labels if you want to use measure %s\n',measure);
  end
end


models = trainedClassifier.models;
classifierType       = trainedClassifier.classifier;
classifierParameters = trainedClassifier.classifierParameters;

lastpos           = size(models,1);
trainingSetInfo   = trainedClassifier.trainingSetInfo;
sortedLabelValues = trainingSetInfo.sortedLabelValues; 
nClasses          = length(sortedLabelValues);
nFeatures         = trainingSetInfo.nFeatures;
nExamples         = size(scores,1);

% 1) Produce matrix with label rankings.
% Each row contains, for one example, all the label values,
% ranked by their decreasing scores. The most probable/highest
% scoring label is thus that on the first column.

  switch classifierType
    
   case {'nbayes','nbayesPooled','nbayes-unitvariance','lda','qda','ldaSVD',...
           'ldaCV','qdaSVD','qdaCV','knn','svmlight','pairwise','neural',...
           'svm','nnets','nbayesPooledResampling','logisticRegression',...
           'logisticRegression2','SMLR','logisticRegressionMulti',}
    trainingSetInfo = models{nClasses+2};
    classPriors     = trainingSetInfo.classPriors;
    logClassConst   = zeros(size(classPriors));  
    
    labelRankings = zeros(nExamples,nClasses);


    % first, obtain a ranking of the scores across all classes (colunns)
    % for each example (row) ( sort(-1*x) sorts x decreasingly)
    [sscores,sindices] = sort( -1*scores, 2 );
    
    % then transform the ranking into a ranking of labels      
    labelRankings = reshape( sortedLabelValues(sindices(:)) ,size(sscores));
    
   otherwise
    fprintf('summarizePredictions: error: classifier %s is not supported\n',classifierType);
    return;
  end
  
  % 2) And now that we have the label ranking matrix, pull out predicted
  % labels for each example from that  
  predictedLabels = labelRankings(:,1);

  % WARNING: this won't work for experiment
  % other than 1ofN if one of the labels is 0
  if isempty(find(sortedLabelValues==0))
    % map numeric label values to indices in the
    % cell array that contains all label values
    for l=1:length(sortedLabelValues)
      label = sortedLabelValues(l);
      labelValuesToLabelIndices{label} = l;
    end
  else
    fprintf('summarizePredictions: WARNING: before you use averageRank makes sure no label values are 0\n'); return;
  end
  
  % Compute total score, output trace with total rankings per class
  nUniqueLabels = length(sortedLabelValues);
  
  % 3) Now switch code depending on the scoring measure desired
  % The code here comes from analysis_runClassifier_dir, and
  % there are other measures there. For the time being, we have
  % accuracy, error and averageRank.

  switch measure
    
   case {'averageRank'}

    % For a given example, measure is the position of the correct
    % class in a ranking of all the label values by their score.
    %
    % 1 - if the correct label has the highest score
    % ...
    % <#classes> - if the correct label has the lowest score
    %
    % The columns in the score array are in the order of sorted
    % label values.
    
    total = 0; % final result - avg rank of the correct label
    trace = cell(2,1); % will contain trace data structures
        

    % counts how many times each label was ranked as another one
    % label 1  <# ranked 1> ... <# ranked k>
    % ...
    % label k  <# ranked 1> ... <# ranked k>

    % counts of how often examples with true label <row> got classified as <col>
    labelCounts  = zeros(nUniqueLabels,nUniqueLabels);
    labelScores  = cell(nExamples,1); % keeps scores of labels
    measureEx    = zeros(nExamples,1); % measure value for example
    
    for t=1:nExamples
      correct  = testLabels(t);
      learnt   = predictedLabels(t);
      % find the position of the true and learnt labels in the ranking of labels
      % by value
      measureEx(t) = find( labelRankings(t,:) == correct );
      correctIdx   = labelValuesToLabelIndices{correct}; % pos of correct
      posIdx       = labelValuesToLabelIndices{learnt};  % pos of learnt
      labelCounts(correctIdx,posIdx) = labelCounts(correctIdx,posIdx) + 1;    
%      fprintf('ex %d\t%d %d\trank=%d\t',t,correct,learnt,pos);
%      fprintf('%d ',labelRankings{t}); fprintf('\n');
    end
    
    % Compute averageRank and normalize to [0,1]
    measureTotal = sum(measureEx) / nExamples; % if all correct, 1
    measureTotal = (measureTotal-1)/(nClasses-1);

    % pack the results and traces
    result{1} = measureTotal;
    result{2} = measureEx;
    trace{1}  = labelRankings;
    trace{2}  = labelCounts;

    
    
   case{'error','accuracy'}

    % Error and accuracy
    % testLabels are the labels of the test examples
    % predictedLabels are those guessed by the classifier    
    nUniqueLabels = length(sortedLabelValues);
    trace  = cell(2,1);
    measureEx = (testLabels == predictedLabels);    
    nright = sum( measureEx );
    ntotal = length( measureEx );
    total  = nright/ntotal;

    % compute labelCounts, the confusion matrix
    % WARNING: this won't work for experiment
    % other than 1ofN if one of the labels is 0
    labelCounts  = zeros(nUniqueLabels,nUniqueLabels);
    
    if isempty(find(sortedLabelValues==0))
      % map numeric label values to indices in the
      % cell array that contains all label values
      for l=1:length(sortedLabelValues)
	label = sortedLabelValues(l);
	labelValuesToLabelIndices{label} = l;
      end
    else
      fprintf('summarizePredictions: WARNING: before you use averageRank makes sure no label values are 0\n'); return;
    end

    for t=1:nExamples
      correct      = testLabels(t);
      learnt       = predictedLabels(t);
      correctIdx   = labelValuesToLabelIndices{correct}; % pos of correct
      posIdx       = labelValuesToLabelIndices{learnt};  % pos of learnt
      labelCounts(correctIdx,posIdx) = labelCounts(correctIdx,posIdx) + 1;    
    end

    % pack the results and traces
    switch measure
     case{'accuracy'}
      result{1} = total;
      result{2} = measureEx;
     case{'error'}
      result{1} = 1 - total;
      result{2} = ~measureEx;
    end
      
    trace{1}  = labelRankings;
    trace{2}  = labelCounts;

   otherwise
    fprintf('summarizePredictions: error: measure %s is not implemented\n',measure);    
  end
  
