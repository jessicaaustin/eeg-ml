% Take a script file and do the experiment described in it
%
% This involves several stages, all preceded by a title on a large
% text box (see "Deal with parameters and the wrapper file" below).
%
%
%
% In:
% - a script file
%
% Out:
% - a directory containing the files used in and the products of
% the experiment such as the classifiers learnt and the measures of
% performance (depends on the classifier and measure picked).
%
% Other functions defined here:
% - idmTransform - creates an IDM of the desired dataType
% - idmToExamples - manages transformation of an IDM into examples for a given experiment
%
% Dependencies:
% - runExperiment depends on MANY functions...
%
% History
% - 23 Aug 02 - fp - created
% 
%
%

function runScript( varargin )  

  present = pwd;  
  l = length(varargin);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %
  % Deal with parameters and the wrapper file
  %  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  ofd = 1;

  switch l
   case{1}
    % passing it a command file, read and set variables from there    
    configFile  = varargin{1}
    subjects    = {};
    rois        = {};
    dataType    = {};
    dataTypeParameters   = {};
    experiment  = {};
    experimentParameters = {};
    classifiers = {};
    datasets    = {};
    group       = {};
    idmPreProcessSteps      = {};
    idmPostProcessSteps     = {};
    examplePreProcessSteps  = {};
    exampleFoldProcessSteps = {};
    doSetup     = 1;
    doTrain     = 1;
    doEval      = 1;
    trainingType = 'L1O'; % by default, leave-1-out
    trainingTypeParameters = {};
    
    ifd = fopen(configFile,'r');
    while 1
      tline = fgetl(ifd);
      if ~ischar(tline), break, end

      try eval(tline);,catch fprintf(1,'runExperiment: error evaluating %s\n',tline); end
    end
    fclose(ifd);

    % Check that file has statements for a number of variables
    % - dataset
    % - experiment
    % - classifiers
    % - dataType
    % - groups
    % - rois
    res = isempty(datasets)+isempty(classifiers)+isempty(experiment)+isempty(dataType)+isempty(group);
    if res > 0
      fprintf('runExperiment: error: please make sure that the wrapper file has the following fields set:\n');
      fprintf('- datasets\n');
      fprintf('- experiment\n');
      fprintf('- classifiers\n');
      fprintf('- dataType\n');
    end
    doConfigFromFile = 1;
    
   otherwise
    fprintf('syntax: runExperiment(<wrapper file>)\nSee template.wrapper for more details\n');
    return;    
  end
  
  nclassifiers = length(classifiers);
  ndatasets    = length(datasets);
  
  % now the parts that vary according to dataset 
  
  for nd=1:1:ndatasets   
    information = mri_reference(datasets{nd});

    % if code does not provide subjects, use mri_reference ones
    if ~doConfigFromFile | ~length(subjects)
      datasetSubjects{nd} = information.subjects;
      datasetMasks{nd}    = information.mask;
    else
      % code provides subjects
      % HACK - assume it is not a dual dataset experiment
      % WARNING: FIX THIS
      datasetSubjects{nd} = subjects;
      datasetMasks{nd}    = mask;      
    end    
  end
      
  % If the wrapper file does not specify which subjects/rois to use
  % we default to the all the ones provided by mri_reference.

  % in the case of joint datasets we assume the first one
  % has all the relevant stuff  
  information = mri_reference(datasets{1}); % {1} IS NOT A BUG
  dataplace   = information.dataplace;
  
  if ~doConfigFromFile | ~length(rois)
    % if file does not provide ROIs
    rois        = information.rois;    
  end
  
  if ~doConfigFromFile | ~length(subjects)
    % if file does not provide subjects
    subjects = information.subjects;
    mask     = information.mask;
  end

  if length(idmPreProcessSteps); idmPreProcess = 1; else; idmPreProcess = 0; end
  if length(idmPostProcessSteps); idmPostProcess = 1; else; idmPostProcess = 0; end
  if length(examplePreProcessSteps); examplePreProcess = 1; else; examplePreProcess = 0; end
  if length(exampleFoldProcessSteps); exampleFoldProcess = 1; else; exampleFoldProcess = 0; end
  dimensions  = information.dimensions;
  dimx        = dimensions(1); dimy = dimensions(2); dimz = dimensions(3);
  nconds      = information.nconds;
    
  % process data to examples figures out datalocs as needed, so
  % this is obsolete. Nevertheless, leave for the time being until
  % fp sorts out the places that depend on this
  dataloc = sprintf('%s/%s',dataplace,datasets{1});

  % assume there is a 1-1 correspondence between the subjects
  % vector of all the datasets in dataset. Guide the loop based
  % on the first one
  nsubjs    = size(subjects,2);
  subjects  = sort(subjects);  
  rois      = sort(rois);
  nrois     = length(rois);

  %
  % A few of the experiment types are special and bypass train/eval
  % Set flags for those.
  
  switch experiment
   case {'metaRandomize'}
    doMeta = 1;
   otherwise
    doMeta = 0;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %
  % Set up data for each experiment
  %  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
	
  % for the time being, only one dataset at once
  dataset = datasets{1};
    
  if doSetup
    for s=1:1:nsubjs
      if mask(s)
      
	subject = subjects{s};
	fprintf(ofd,'setting up data for subject %s\n',subject);
	
	msg = sprintf('\tcreating datasets for experiment %s\n',experiment);	  
	fprintf(ofd,msg);
	
	if isequal(group,'separate') useSeparateROIs=1; else useSeparateROIs=0; end

	% 0) Create directory to keep stuff
	% if experiment names can be parameterized need more information
	edir = nameDataDir(experiment,subject,group);
	s=mkdir(edir);
	cd(edir);

	if ~doMeta
	  % This is the normal case - metaExperiments bypass this
	  
	  % 1) Retrieve IDMs for all the rois in this subject
	  %    They will come already merged into a single large ROI
	  
	  [is,ds,ms] = loadSubjectdataMult(dataset,subject,rois,useSeparateROIs);
	  
	  % 2) Create examples and save into files the next stage can use	
	  % Create examples from each IDM
	  % If useSeparateROIs is 0 we have a single ROI here
	  
	  for r=1:length(ms)
	    roi = ms{r}.roi;
	    
	    % a) apply IDM space pre processing
	    if idmPreProcess
	      [is{r},ds{r},ms{r}] = transformIDM_preProcess(is{r},ds{r},ms{r},idmPreProcessSteps);
	    end
	    
	    % b) apply datatype related IDM transformations	  
	    % This part depends on what study the data came from and
	    % the data type selected
	    [is{r},ds{r},ms{r}] = transformIDM_dataType(is{r},ds{r},ms{r},dataType,dataTypeParameters,experiment,experimentParameters,datasets);

	    % c) apply IDM space post processing
	    if idmPostProcess
	      [is{r},ds{r},ms{r}] = transformIDM_preProcess(is{r},ds{r},ms{r},idmPostProcessSteps);
	    end
	    	    
	    % d) IDM -> examples (depends on the experiment)
	    [examples,labels,expInfo] = idmToExamples(is{r},ds{r},ms{r},dataType,experiment,datasets);	  
	    % e) TODO: apply example space preprocessing
	    if examplePreProcess
	      [examples,labels] = transformExamples_preProcess(examples,labels,examplePreProcessSteps);
	    end
	    	    
	    % f) now output the following into a file:	  
	    %   examples - contains examples
	    %   labels   - contains labels
	    %   dataset+dataType+subject+roi+experiment
	    %   examplePreProcessInfo
	    %   dataPreProcessInfo
	    %   expInfo  - contains information about the
	    %   IDM->examples process
	    expInfo.experiment = experiment;
	    expInfo.meta       = ms{r};
	    
	    exampleFile = sprintf('examples.%s.%s.%s.mat',experiment,roi,dataType);
	    
	    cmd = sprintf('save %s examples dataset expInfo labels experiment dataType subject roi useSeparateROIs idmPreProcessSteps idmPostProcessSteps examplePreProcessSteps',exampleFile);
	    eval(cmd);

	    % if what we have is a metaExperiment, just go ahead
            % and run the entire thing from here - a bit of a hack...
	    
	    % pack the relevant parts of our split
	    % subjects    = {};
	    % rois        = {};
	    % dataType    = {};
	    % experiment  = {};
	    % classifiers = {};
	    % datasets    = {};
	    % group       = {};
	    % idmPreProcessSteps     = {};
	    % idmPostProcessSteps    = {};
	    % examplePreProcessSteps = {};
	    % dataTypeParameters   = {};
	    % experimentParameters = {};    
	    % doSetup     = 1;
	    % doTrain     = 1;
	    % doEval      = 1;	    

	  end; % for ROIS
	end;% if doMeta
	  
	cd ..;

      end ;% if subject is in the mask
    end ;% for over all the subjects
  end

  % at this point, we have the IDMs for one or more ROIs in is,ds,ms
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %
  % Train classifiers and produce scores
  %  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  if doTrain & ~doMeta

    fprintf(ofd,'\n\n-----------------------------------\n\n');
      
    for s=1:1:nsubjs
      if mask(s)
	
	subject = subjects{s};
	fprintf(ofd,'Training classifiers for subject %s using %s\n',subject,trainingType);
    
	edir = nameDataDir(experiment,subject,group);
	cd(edir);

	% Train classifier for each set of examples

	%	for r=1:length(rois)
	%	  roi = rois{r};
	% temporary fix to allow training without setting
	
	dummy = dir(sprintf('examples.%s.*.%s.mat',experiment,dataType)); 
	rest  = dummy(1).name;
	[prefix,rest] = strtok(rest,'.');
	[prefix,rest] = strtok(rest,'.');
	[roi,rest] = strtok(rest,'.');
	
	for cl=1:1:nclassifiers
	  classifier = classifiers{cl};
	  
	  fprintf(ofd,'\tclassifier %s\n',classifier);
	  
	  % load examples
	  exampleFile = sprintf('examples.%s.%s.%s.mat',experiment,roi,dataType);
	  load(exampleFile); 
	  sortedLabelValues = sort(unique(labels));
	  nClasses    = length(sortedLabelValues); 
	  
	  % we don't figure example/feature things yet, as we may
	  % have preprocessing to apply
	  
	  switch trainingType
	    % trainingType names are provisional
	    
	   case {'halfAndHalf'}
	    % split the data in half - (using the fact that the number is even)
	    nExamples = size(examples,1);
	    nFeatures = size(examples,2);
	    ntotal    = size(examples,1);
	    oindices  = 1:2:(ntotal-1);
	    eindices  = 2:2:ntotal;
	    trainExamples = examples(oindices,:);
	    trainLabels   = labels(oindices,1);
	    testExamples  = examples(eindices,:);
	    testLabels    = labels(eindices,1);
	    % train a classifier
	    [models] = trainClassifier(trainExamples,trainLabels,'nbayes');  
	    % apply a classifier
	    [scores] = applyClassifier(testExamples,models,'nbayes');
	    
	   case {'givenTrainTest'}
	    % expInfo knows which examples are in the training
	    % set and which are in the testing set
	    nExamples = size(examples,1);
	    nFeatures = size(examples,2);
	    scores    = zeros(nExamples,nClasses);
	    
	    % 1st block of a condition as training set
	    % 2nd block of a condition as testing set
	    trainingMask  = expInfo.training;
	    trainIdx      = find(trainingMask);
	    testIdx       = find(~trainingMask);
	    trainExamples = examples(trainIdx,:);
	    trainLabels   = labels(trainIdx,:);
	    testExamples  = examples(testIdx,:);
	    testLabels    = labels(testIdx,:);	      
	    % train a classifier
	    [models1] = trainClassifier(trainExamples,trainLabels,'nbayes');  
	    % apply a classifier
	    [scores1] = applyClassifier(testExamples,models1,'nbayes');
	    scores(testIdx,:) = scores1;
	    
	    % 1st block of a condition as testing set
	    % 2nd block of a condition as training set
	    trainingMask = ~expInfo.training;
	    trainIdx     = find(trainingMask);
	    testIdx      = find(~trainingMask);
	    trainExamples = examples(trainIdx,:);
	    trainLabels   = labels(trainIdx,:);
	    testExamples  = examples(testIdx,:);
	    testLabels    = labels(testIdx,:);	      
	    % train a classifier
	    [models2] = trainClassifier(trainExamples,trainLabels,'nbayes');  
	    % apply a classifier
	    [scores2] = applyClassifier(testExamples,models2,'nbayes');
	    scores(testIdx,:) = scores2; 
	    
	    % finally compute a model over all the data
	    [models] = trainClassifier(examples,labels,'nbayes');
	    
	   case {'kFoldCV'}
	    % kFold cross validation
	    nExamples = size(examples,1);
	    nFeatures = size(examples,2);
	    nFolds    = trainingTypeParameters{1};
	    
	    fprintf('doing %d-fold cross validation\n',nFolds);
	    
	    % Shuffle examples prior to separating them into folds
	    % This is a quick "solution" to the problem of ending
            % up with fold training sets that don't have examples
            % of a given class. We could try a more equitative distribution
	    ensemble    = sortrows([(1:nExamples)',(randperm(nExamples))',labels,examples],2);
	    newExamples = ensemble(:,4:nFeatures+3);
	    newLabels   = ensemble(:,3);
	    newOrder    = ensemble(:,1);

	    if exampleFoldProcess
	      [models,newScores] = trainClassifier_kFoldCV(newExamples,newLabels,expInfo,classifier,dataType,nFolds,exampleFoldProcessSteps);
	    else
	      [models,newScores] = trainClassifier_kFoldCV(newExamples,newLabels,expInfo,classifier,dataType,nFolds);
	    end
	    
	    % sort scores back into the original order
	    ensemble = sortrows([newOrder,newScores],1);
	    scores   = ensemble(:,2:nClasses+1);
	    
	    clear newExamples newScores newLabels ensemble
	    
	   otherwise
	    % by default, do leave-1-out
	    nExamples = size(examples,1);
	    nFeatures = size(examples,2);
	    
	    if exampleFoldProcess
	      [models,scores] = trainClassifier_L1O(examples,labels,expInfo,classifier,dataType,exampleFoldProcessSteps);
	    else
	      [models,scores] = trainClassifier_L1O(examples,labels,expInfo,classifier,dataType);
	    end
	  end
	    
	  % save learnt models and scores
	  modelFile = sprintf('models.%s.%s.%s.%s.mat',experiment,roi,dataType,classifier);
	  scoreFile = sprintf('scores.%s.%s.%s.%s.mat',experiment,roi,dataType,classifier);
	  
	  cmd = sprintf('save %s models dataset experiment dataType subject roi useSeparateROIs',modelFile); eval(cmd);
	  cmd = sprintf('save %s scores dataset experiment dataType subject roi useSeparateROIs',scoreFile); eval(cmd);
	  
	end ;% for classifiers	  
%	end; % for rois
cd ..;
 
      end ;% if mask
    end ;% for subjects
    
  end; % if doRun
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %
  % Evaluate results
  %  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  if doEval & ~doMeta
    
    fprintf(ofd,'\n\n-----------------------------------\n\n');
    
    for s=1:1:nsubjs
      if mask(s)
	
	subject = subjects{s};
	fprintf(ofd,'Evaluating subject %s\n',subject);
    
	edir = nameDataDir(experiment,subject,group);
	cd(edir);

	% Train classifier for each set of examples
	%	for r=1:length(rois)
	%	  roi = rois{r};    
	  
	  % temporary fix to allow evaluating without setting up
	  
	  dummy = dir(sprintf('examples.%s.*.%s.mat',experiment,dataType)); 
	  rest  = dummy(1).name;
	  [prefix,rest] = strtok(rest,'.');
	  [prefix,rest] = strtok(rest,'.');
	  [roi,rest] = strtok(rest,'.');
	  
	  for cl=1:1:nclassifiers
	    classifier = classifiers{cl};
	  
	    fprintf(ofd,'\tclassifier %s\n',classifier);

	    % load models and scores
	    exampleFile = sprintf('examples.%s.%s.%s.mat',experiment,roi,dataType);
	    modelFile = sprintf('models.%s.%s.%s.%s.mat',experiment,roi,dataType,classifier);
	    scoreFile = sprintf('scores.%s.%s.%s.%s.mat',experiment,roi,dataType,classifier);
	    load(exampleFile);
	    load(modelFile); 
	    load(scoreFile); 

	    switch experiment
	     case {'metaExperiment','metaContrastConds','metaOneOther','metaCompareBlocks'}
	      % TODO
	     otherwise
	      [result,predictedLabels,confusionMatrix] = summarizePredictions(examples,models,scores,classifier,'averageRank',labels);
	    end
	    
	    % output information into trace files
	    % TODO: move this out into a new function
	    traceFile = sprintf('trace.%s.%s.%s.%s.txt',experiment,roi,dataType,classifier);
	    outfd   = fopen(traceFile,'w');
	    nExamples = size(confusionMatrix{1},1);
	    nClasses  = size(confusionMatrix{1},2);
	    scores = double(scores);
	    
	    fprintf(outfd,'%% Average Rank\n');
	    fprintf(outfd,'%% %1.3f\n',result{1});	    
	    fprintf(outfd,'%% ');
	    fprintf(outfd,'Pred\tTrue\t');
	    fprintf(outfd,'Ranks\t');
	    for k=2:nClasses; fprintf(outfd,'- '); end
	    fprintf(outfd,'Scores\t');
	    for k=2:nClasses; fprintf(outfd,'- '); end
	    fprintf(outfd,'\n');
	    
	    for n=1:nExamples
	      fprintf(outfd,'%d\t',predictedLabels(n,1));
	      fprintf(outfd,'%d\t',labels(n,1));    
	      fprintf(outfd,'%d ',confusionMatrix{1}(n,:));
	      fprintf(outfd,'\t');
	      fprintf(outfd,'%1.1f ',scores(n,:));
	      fprintf(outfd,'\n');
	    end  
	    fclose(outfd);    

	    % output results
	    resFile = sprintf('result.%s.%s.%s.%s.txt',experiment,roi,dataType,classifier);
	    outfd   = fopen(resFile,'w');
	    fprintf(outfd,'%1.3f',result{1});
	    fclose(outfd);
	    
	    % output confusion matrix
	    cmatFile = sprintf('cmatrix.%s.%s.%s.%s.txt',experiment,roi,dataType,classifier);
	    outfd   = fopen(cmatFile,'w');
	    nExamples = size(confusionMatrix{1},1);
	    nClasses  = size(confusionMatrix{1},2);
	    
	    for n=1:nClasses
	      fprintf(outfd,'%d\t',confusionMatrix{2}(n,:));
	      fprintf(outfd,'\n');
	    end  
	    fclose(outfd);
	    
	  end ;% for classifiers	  
	%end; % for rois
	
	cd ..;
	
      end ;% if mask
    end ;% for subjects
  end; % if eval
    
  
  
%
% Creates an IDM with the desired dataType (each dataType specifies
% a pipeline of transformations)
%
% dataTypeParameters and experimentParameters are cell arrays of
% parameters that a few of the dataTypes require (e.g. how many
% voxels to keep when cropping to active voxels in the "active"
% dataType).
%

function [info,data,meta] = transformIDM_dataType(info,data,meta,dataType,dataTypeParameters,experiment,experimentParameters,datasets)

  % First do transformations that will determine how many trials
  % and voxels there will ultimately be (separate conditions,
  % average trials, pick active voxels, etc)
  
  switch dataType
    
   case{'full'}
    
   case {'unfolded'}
    % Turn each time point in a trial into a new trial with length
    % 1, same condition as the original trial
    [info,data,meta] = transformIDM_unfold(info,data,meta);

   case {'active'}
    % Turn into IDM with the union of the <nToKeep> most active voxels in each condition
    % Get <nToKeep> from the dataType parameters
    nToKeep = dataTypeParameters{1};    
    [info,data,meta] = transformIDM_selectActiveVoxels(info,data,meta,nToKeep);
    
   case{'unfoldedSeparateBlocks'}
    % Turn trials with the same condition into separate new conditions
    [info,data,meta] = transformIDM_separateBlocks(info,data,meta);
    % Turn each time point in a trial into a new trial with length
    % 1, same condition as the original trial
    [info,data,meta] = transformIDM_unfold(info,data,meta);
    
   case{'unfoldedWindowAvg'}
    % average all the trials in the same condition into a single one
    [info,data,meta] = transformIDM_pairwiseAvg(info,data,meta);
    % Turn each time point in a trial into a new trial with length
    % 1, same condition as the original trial
    [info,data,meta] = transformIDM_unfold(info,data,meta);
    
   case{'unfoldedAvgTrial'}
    % average all the trials in the same condition into a single one
    [info,data,meta,avgstd] = transformIDM_averageOverTrials(info,data,meta);
    % Turn each time point in a trial into a new trial with length
    % 1, same condition as the original trial
    [info,data,meta] = transformIDM_unfold(info,data,meta);
    
   case{'unfoldedActive','unfoldedActiveSeparateBlocks'}   
    if isequal(dataType,'unfoldedActiveSeparateBlocks')
      % Turn trials with the same condition into separate new conditions
      [info,data,meta] = transformIDM_separateBlocks(info,data,meta);
    end
    % This will be an IDM with the union of the <nToKeep> most
    % active voxels in each condition in each ROI
    % Get <nToKeep> from the dataType parameters
    nToKeep = dataTypeParameters{1};
    
    [info,data,meta] = transformIDM_selectActiveVoxels(info,data,meta,nToKeep);
    % Turn each time point in a trial into a new trial with length
    % 1, same condition as the original trial
    [info,data,meta] = transformIDM_unfold(info,data,meta);

   case{'flat'}
    % now flatten each trial to a single time point
    [info,data,meta] = mri_flattenTrialFilter(info,data,meta);
    
   case{'snapshot'}
    [info,data,meta] = mri_snapshot(tSlice,info,data,meta);
    
   otherwise
    % this is a distilled data abstraction
    % which has been pre computed
    fprintf('runExperiment: error: dataType %s is not supported\n',dataType);
    return;
  end

  % gather information about number,type and length of trials.
  % also figure out how much of a trial to use: [trialBegin,trialEnd]
  dataset = datasets{1};
  
  [ntrials,nvoxels,nconds,minTrialLenCond,ntrialsCond,trialBegin,trialEnd] = mri_infoTrials( info,data,meta,dataset);
  
  information = mri_reference(dataset);
  dimensions  = information.dimensions;
  dimx        = dimensions(1); dimy = dimensions(2); dimz = dimensions(3);
  
  
%
% Transform an IDM into a set of examples for the desired experiment
%

function  [examples,labels,expInfo] = idmToExamples(info,data,meta,dataType,experiment,datasets)

  dataset = datasets{1};
  
  switch dataset
    
    %%%%%%%%%% STARPLUS EXPERIMENTS %%%%%%%%%%%%
   case {'data-starplus-sp','data-starplus-ps'}
    
    switch experiment
     otherwise
      fprintf(1,'\terror: experiment %s.%s is not available\n',dataset,experiment);return;
    end

    %%%%%%%%%% CATEGORIES EXPERIMENTS %%%%%%%%%%%%
   case {'data-categories'}
    
    switch experiment
     case {'fixation'}
      [examples,labels,expInfo] = idmToExamples_fixation(info,data,meta,dataType);
     case {'1ofN'}
      [examples,labels,expInfo] = idmToExamples_1ofN(info,data,meta,dataType);
     case {'1ofNblockFromBlock'}
      [examples,labels,expInfo] = idmToExamples_1ofNbFb(info,data,meta,dataType);      
      
     otherwise
      fprintf(1,'\terror: experiment %s.%s is not available\n',dataset,experiment);return;
    end
    
    %%%%%%%%%% SYNTACTIC AMBIGUITY EXPERIMENTS %%%%%%%%%%%%
   case {'data-syntamb2'}
    
    switch experiment
     otherwise
      fprintf(1,'\terror: experiment %s.%s is not available\n',dataset{1},experiment);
    end
    
   otherwise
    fprintf(1,'\terror: dataset %s is not supported\n',dataset);return;
  end

  expInfo.meta = meta;
  

  
% Apply preprocessing steps to IDM  
%
% The wrapper file contains a cell array specifying a sequence of
% transformations (one per cell). Each transformation comprises:
% 1) function to apply to IDM
% 2) arguments required in addition to IDM
%
% This function goes through the cell array and applies each transformation

function [info,data,meta] = transformIDM_preProcess(info,data,meta,processSteps)
    
  % go through all the steps  
  nSteps = length(processSteps);
  
  processSteps
  
  for s=1:nSteps

    step = processSteps{s};
    stepName = step{1};
    stepArgs = step{2};
    nArgs    = length(stepArgs);
    fprintf('IDM preprocess step %d - %s\n',s,stepName);
    
    newArgs  = cell(nArgs+3,1);
    newArgs{1} = info; newArgs{2} = data; newArgs{3} = meta;
    for a=1:nArgs; newArgs{a+3} = stepArgs{a}; end
    
    cmd = sprintf('[info,data,meta]=%s(newArgs);',stepName);
    eval(cmd);
    
    if 0
      % uncomment to have printouts of what is being called in each
      % step, and with what arguments      
      for a=1:nArgs
	if isnumeric(stepArgs{a})
	  isfloat = stepArgs{a} - round(stepArgs{a});
	  if isfloat
	    fprintf('\t%f\n',stepArgs{a});
	  else
	    fprintf('\t%d\n',stepArgs{a});
	  end	
	elseif ischar(stepArgs{a})
	  fprintf('\t%s\n',stepArgs{a});
	else
	  fprintf('\tCannot handle type of argument %d\n',a);
	end
      end
    end
  
  end

  
% Apply preprocessing steps to an array of examples
%
% The wrapper file contains a cell array specifying a sequence of
% transformations (one per cell). Each transformation comprises:
% 1) function to apply to examples+labels
% 2) arguments required in addition to examples+labels
%
% This function goes through the cell array and applies each transformation

function [examples,labels] = transformExamples_preProcess(examples,labels,processSteps)
    
% go through all the steps  
  nSteps = length(processSteps);
  
  processSteps
  
  for s=1:nSteps

    step = processSteps{s};
    stepName = step{1};
    stepArgs = step{2};
    nArgs    = length(stepArgs);
    fprintf('Example preprocess step %d - %s\n',s,stepName);
    
    newArgs  = cell(nArgs+2,1);
    newArgs{1} = examples; newArgs{2} = labels;
    for a=1:nArgs; newArgs{a+2} = stepArgs{a}; end
    
    cmd = sprintf('[examples,labels]=%s(newArgs);',stepName);
    eval(cmd);
    
    if 0
      % uncomment to have printouts of what is being called in each
      % step, and with what arguments      
      for a=1:nArgs
	if isnumeric(stepArgs{a})
	  isfloat = stepArgs{a} - round(stepArgs{a});
	  if isfloat
	    fprintf('\t%f\n',stepArgs{a});
	  else
	    fprintf('\t%d\n',stepArgs{a});
	  end	
	elseif ischar(stepArgs{a})
	  fprintf('\t%s\n',stepArgs{a});
	else
	  fprintf('\tCannot handle type of argument %d\n',a);
	end
      end
    end
  
  end

% used wherever a data directory name is needed - this keeps things
% coherent

function [edir] = nameDataDir(experiment,subject,group)
  
% Create a directory name
% some experiments can be parameterized and hence get different
% names for every setting
  switch experiment
   case {'contrastPairConds'}
    % not yet
   case {'contrastOneOther'}
    % not yet
   otherwise
    experimentName = experiment;
  end

  edir = sprintf('datafor_%s_%s_%s',experimentName,subject,group);

  
  
%
% A few simple IDM transforms not worth putting into a separate function
% (at least yet).

% Randomizes the block label assignments

function [info,data,meta] = transformIDM_randomizeBlocks(notvarargin)

  varargin = notvarargin;
  l = length(varargin);
  info = varargin{1}; data = varargin{2}; meta = varargin{3};
  if l > 3; seed = varargin{4}; else seed = 4791; end
  
  % gather information about number,type and length of trials and other dataset parameters
  [ntrials,nvoxels,nconds,minTrialLenCond,ntrialsCond,trialBegin,trialEnd] = mri_infoTrials( info,data,meta,meta.study );

  % establish a new randomized labeling for the blocks  
  rand('seed',seed);
  blockList = [];
  blockIdx  = 1;
  
  for nt=1:ntrials
    if info(nt).cond > 1
      blockList(blockIdx) = info(nt).cond;
      blockIdx = blockIdx + 1;    
    end
  end
  
  nblocks  = blockIdx - 1;  
  ensemble = sortrows(([randperm(nblocks);blockList])',1);
  newBlockList = ensemble(:,2)';
  
  disp(blockList);
  disp(newBlockList);
  
  % now relabel  
  blockIdx = 1;
  for nt=1:ntrials
    if info(nt).cond > 1
      info(nt).cond = newBlockList(blockIdx);
      blockIdx = blockIdx + 1;    
    end
  end
  
  
% Also randomizes block label assignments, but:
% - assumes there are two blocks per condition
% - assumes condition numbers <= # classes + 1
% - no two blocks that had the same condition originally will
%   have it in the new assignment

function [info,data,meta] = transformIDM_randomizeBlocksP(notvarargin)

  varargin = notvarargin;
  l = length(varargin);
  info = varargin{1}; data = varargin{2}; meta = varargin{3};
  if l > 3; seed = varargin{4}; else seed = 4791; end
  
  % gather information about number,type and length of trials and other dataset parameters
  [ntrials,nvoxels,nconds,minTrialLenCond,ntrialsCond,trialBegin,trialEnd] = mri_infoTrials( info,data,meta,meta.study );

  % establish a new randomized labeling for the blocks  
  rand('seed',seed);
  blockList = [];
  blockIdx  = 1;  
  
  for nt=1:ntrials
    if info(nt).cond > 1
      blockList(blockIdx) = info(nt).cond;
      blockIdx = blockIdx + 1;    
    end
  end
  
  nblocks  = blockIdx - 1;  
  nClasses = nblocks/2;
  
  if floor(nClasses) ~= nClasses
    fprintf('transformIDM_randomizeBlocksP: ERROR: 2 blocks per condition required\n');
    return;
  end
  
  notGood = 1;
  idx = 1;
  while( notGood )    
    fprintf('\t\ttesting permutation %d\n',idx);
    
    pairCounts = zeros(nClasses,nClasses);
    
    ensemble = sortrows(([randperm(nblocks);blockList])',1);
    newBlockList = ensemble(:,2)';
    disp(blockList);
    disp(newBlockList);
    
    % test the assignment
    indices = sub2ind([nClasses+1 nClasses+1],blockList,newBlockList);
    if length(indices) == length(unique(indices)); notGood = 0; end
    disp(notGood);
    idx = idx + 1;
  end
  
  % now relabel  
  blockIdx = 1;
  for nt=1:ntrials
    if info(nt).cond > 1
      info(nt).cond = newBlockList(blockIdx);
      blockIdx = blockIdx + 1;    
    end
  end

  
% template IDM preprocessing code

function [info,data,meta] = transformIDM_dummy(notvarargin)

  varargin = notvarargin;
  l = length(varargin);
  info = varargin{1}; data = varargin{2}; meta = varargin{3};
  number = varargin{4};
  
  fprintf('transformIDM_dummy with number=%d\n',number);
  
% template Example preprocessing code
  
function [examples,labels] = transformExamples_dummy(notvarargin)

  varargin = notvarargin;
  l = length(varargin);
  examples = varargin{1}; labels = varargin{2};
  number = varargin{3};
  
  fprintf('transformExamples_dummy with number=%d\n',number);