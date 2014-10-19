% Do leave 1 out training and testing
% Output discriminative/generative models trained on the full dataset.
%
% WARNING: it only does leave 1 out right now, but should be easy
%          to modify it to do k-fold cross validation
%
% In:
% - examples, labels, classifier, experiment and dataType the
%   examples came from.
%
% Out:
% - a cell array containing models (more details below in NOTES)
% - scores array for all the examples
%
% Notes:
% In addition to the examples+labels, we need to know which
% experiment and dataType the examples came from. One example
% of why this is needed is the "unfolded" data type, where
% a certain number of examples preceding and following those in
% the test fold in time during the experiment need to be excluded
% from the training set. Again in this circumstance, if the data
% comes from an IDM where each example has been replaced by the
% average of itself with its neighbours then the window of
% exclusion will have to be larger. The circumstance is conveyed
% via the data type.
% 
% The models cell array that is output:
%
% models = cell(nClasses+3,1);
%
% Cells 1 to nclasses - contain the generative model for each class
% Cell nclasses+1 - contains the discriminative model, if any
% Cell nclasses+2 - contains information about training set
% Cell nclasses+3 - additional information, if any
%
% Generative model - each cell contains a cell array
% where each cell has one parameter, e.g.
%
%  models{1}{1} = mean class 1
%  models{1}{2} = stdev class 1
%  models{1}{3} = prior prob class 1
%
% Discriminative model - a cell array of sets of weights
%  models{nclasses+1} = cell(1,1);
%  
% Information
%  models{nclasses+2} = cell(1,1);
% containing
%  trainingSetInfo.nExamples (for the time being)
%
% Dependencies:
% - trainLeave1Out depends on classifierNBayes
%
% History
% - 09 Aug 02 - fp - adapted from fmri_devel/Classification/analyser_runClassifier_l1out
%
%

function [models,scores] = trainClassifierL1O( varargin )

  l = length(varargin);
  if l < 4
    fprintf(1,'syntax: trainClassifier(examples,labels,classifier,dataType)\n');
    return;
  end
  
  examples       = varargin{1};
  labels         = varargin{2};
  classifier = varargin{3};
  dataType       = varargin{4};
  
  % retrieve the rest of the information and set up parameters
  nExamples = size(examples,1);
  nFeatures = size(examples,2);
  nLabels   = size(labels,2);
  sortedLabelValues = sort(unique(labels));
  nClasses    = length(sortedLabelValues);
      
  %
  % Simple leave 1 out loop
  %

  % 1) examples are assigned to folds
  % 2) for each testing fold:
  %   a)  run classifier on the remaining data (minus whatever needed to be excluded)
  %   b)  take the model output by the classifier and use it to
  %   produce scores
  % 
  
  % Assign each example to a fold
  % as this is leave-1-out, fold assignment is 1:1:nExamples
  % anything else would require sorting examples by fold assignment
  nFolds    = nExamples;
  foldsize  = nExamples/nFolds;
  results   = zeros(nFolds,1);
  
  foldAssignment = 1:1:nExamples;
  featureCounts  = zeros(nFeatures,1);
  
  nTest  = foldsize;
  nTrain = (nFolds-1)*foldsize;
  
  % run the actual loop
  
  % initialize data structures (so that the first clear does not
  % complain)
  testSet     = zeros(nTest,nFeatures);
  testLabels  = zeros(nTest,nLabels);
  trainSet    = zeros(nTrain,nFeatures);
  trainLabels = zeros(nTrain,nLabels);
  
  foldModels   = cell(nFolds,1);% will contain the models produced for each fold ran
  learntLabels = zeros(nExamples,nLabels); % will contain learnt labels
  scores       = ones(nTest,nClasses)/nClasses; % scores for each class
  foldPreProcessInfo = cell(nFolds,1); % preprocessing info
    

  for k = 1:1:nFolds

    fprintf('\tdoing fold %d\n',k);
    %
    % Set up examples into train+test sets
    % 
    clear testSet testLabels trainSet trainLabels;
    
    testSet     = zeros(nTest,nFeatures);
    testLabels  = zeros(nTest,nLabels);
    
    % copy this fold into the test set
    idxB = 1 + (k-1)*foldsize; % first pos in fold
    idxE = k*foldsize;         % last pos in fold
    testSet    = examples(idxB:1:idxE,:);
    testLabels = labels(idxB:1:idxE,:);

    % figure out how many training examples we will get in this case    
    idxEprev = idxB-1; % last position in previous fold (or <= 0)
    idxBnext = idxE+1; % first position in the next fold
    continueAt = 1;

    % If required by the dataType, eliminate a window of examples around
    % current fold (see Notes at the beginning of the code for
    % details of why and when this happens)
    
    switch dataType
     
     case {'unfolded','unfoldedWindowAvg','unfoldedSeparateBlocks','unfoldedAvgTrial','unfoldedActive','unfoldedActiveSeparateBlocks'}

      switch dataType
       case {'unfolded','unfoldedSeparateBlocks','unfoldedAvgTrial','unfoldedActive','unfoldedActiveSeparateBlocks'}
	rmwindow = 6;
       otherwise
	rmwindow = 8;
      end
      
      idxEprev = idxEprev - rmwindow; 
      idxBnext = idxBnext + rmwindow;
     otherwise
      % some 
    end

    % at this point
    % idxEprev - last example before fold to go into the training set
    % idxBnext - first example after fold to go into the training set

    if idxEprev > 0; nbefore = idxEprev; else; nbefore = 0; end;
    if idxBnext > nExamples; nafter = 0; else; nafter = nExamples-idxBnext+1;end
    
    nTrain      = nbefore + nafter;
    trainSet    = zeros(nTrain,nFeatures);
    trainLabels = zeros(nTrain,nLabels);
    
    % copy previous folds into the training set
    if idxEprev > 0
      trainSet(1:1:idxEprev,:)    = examples(1:1:idxEprev,:);
      trainLabels(1:1:idxEprev,:) = labels(1:1:idxEprev,:);
      continueAt = idxEprev + 1;
    end
    
    % copy subsequent folds into the training set
    trainSet(continueAt:1:nTrain,:)    = examples(idxBnext:1:nExamples,:);	
    trainLabels(continueAt:1:nTrain,:) = labels(idxBnext:1:nExamples,:);    
          
    %
    % Run classifiers and store models and scores
    %
    % idxB - beginning of this fold, idxE - end of this fold
    %
    % Each classifier call over a fold returns learnt models. Those
    % can then be used to obtain scores over the labels of the
    % examples in the fold
    %

    switch classifier

     case {'nbayes'}
      % train the classifier and produce models
      [models] = trainClassifier(trainSet,trainLabels,'nbayes');
      % use the models to score labels for the test set
      [foldScores] = applyClassifier(testSet,models,'nbayes');
      %	[classifierLabels,models,logClassProbs] =
      %	classifier_nbayes( trainSet, trainLabels, testSet);
      %	% what used to be here	
     otherwise
      fprintf('trainLeave1Out: error: classifier %s is not supported\n',classifier);
    end ;% switch over classifiers
  
    % store scores and models

    scores(idxB:idxE,:) = foldScores;
    foldModels{k}       = models; % store learnt models
  
  end ;% for over folds

  %
  % now run the classifier over all the data available to get the final models
  %
  % This would be the place to do more complicated combinations of the models
  % learnt for each fold.
    
  trainSet = examples;
  testSet  = examples;
    
  switch classifier
    
   case {'nbayes'}
    [models] = classifierBayes(examples,labels,'nbayes');    
   otherwise
    fprintf('trainLeave1Out: error: classifier %s is not supported\n',classifier);
  end