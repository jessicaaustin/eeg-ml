%
% Train a SVM-light classifier using cross validation within
% the training set to decide the value of the C parameter.
% Uses the SVM-MEX and SVM-light packages, stored in the SVM subdirectory
%
% In:
% - examples
% - labels
% - optionally, a classifierParameters cell array, containing:
%   - kernel type - 'linear'|'polynomial'|'RBF'|'sigmoid','userDefined'
%     (default 'linear')
%   - kernel parameters - a vector of values which are options for each kernel: 
%                         []       for a linear kernenl
%                         [d]      for polynomial kernel
%                         [gamma]  for RBF
%                         [s, c]   for tanh kernel,
%                         'string' (instead of a vector) for a user-defined kernel
%     (default is no parameters, as kernel default is 'linear')
%   - C parameter - 'svmlightDefault' (use the svmlight default value)
%                   'crossValidation' (use CV to set it)
%                   <value>           (use a specified value)
%     (default is svmlightDefault)
%   - normalize the data - make each feature mean 0 and stdev 1
%     (default is 0)
%
% Out:
% - a learnt SVM model
%
% Example:
% % use linear kernel, default C parameter setting
% - classifierParameters = {'linear',[]};
% % use quadratic kernel, default C parameter setting
% - classifierParameters = {'polynomial',[2]}; % use default parameter setting
% % use linear kernel, use cross-validation to set C
% - classifierParameters = {'linear',[],'crossValidation'};
%
% Notes:
% - works only on two class problems (train/applyClassifier.m handle multiclass)
%
% History
% - 2005 Sep 11 - fpereira - created
%
  
function [model] = classifierSVMlight( varargin )

DEBUG = 0;
this = 'classifierSVMlight';

%
% Process parameters
%
 
l = length(varargin);
if l < 1; help classifierSVMlight; return; end

examples = varargin{1}; [nExamples,nFeatures] = size(examples);
labels   = varargin{2}; sortedLabelValues = unique(labels); nLabels = length(sortedLabelValues);

classifierParameters = {};
if l > 2
  classifierParameters = varargin{3};
end

clear varargin; % it can be very large

if nLabels > 2
  fprintf('%s: this function only supports 2 class problems,call trainClassifier.m for multiclass problems\n',this);pause;return;
end

% defaults
kernel        = 'linear';
kernelParams  = [];
optimC        = 'svmlightDefault';  % let SVMlight pick the default
normalizeData = 0;

ncp = length(classifierParameters);

% override them
if ncp
  
  kernel       = classifierParameters{1}; % kernel type
  
  if ncp > 1;     kernelParams = classifierParameters{2}; % a parameter string
    if ncp > 2;   optimC = classifierParameters{3};
      if ncp > 3; normalizeData = classifierParameters{4};
      end
    end
  end
end

%% create command string to pass to code

% 1) incorporate kernel type
switch kernel
 case {'linear'}
  kernelString = sprintf('-t 0');    
 case {'polynomial'}
  kernelString = sprintf('-t 1 -d %d',kernelParams(1));    
 case {'RBF'}
  kernelString = sprintf('-t 2 -g %d',kernelParams(1));
 case {'sigmoid'}
  kernelString = sprintf('-t 3 -s %d %d',kernelParams(1),kernelParams(2))';
 case {'userDefined'}
  kernelString = sprintf('-t 4 %s',kernelParams);
 otherwise
  fprintf('%s: kernel %s is not supported\n',this);pause;return
end

% make SVMlight silent
kernelString = sprintf('%s ',kernelString);

% 2) add extra parameters to the argument string, if required
useCV = 0;
switch optimC
 case {'svmlightDefault'}
  % the code will do this by default
 case {'crossValidation'}
  % the cross validation code will take care of putting values in the string
  useCV = 1;
 otherwise
  % an actual value
  if isnumeric(optimC)
    kernelString = sprintf('%s -c %s',kernelString,num2str(optimc));
  else
    fprintf('%s: optimC must be a numeric value\n',this);pause;return
  end
end

% transform the labels into +1 and -1
mask = (labels == sortedLabelValues(1));
indices = find(mask);
labels(indices) = 1;
indices = find(~mask);
labels(indices) = -1;

%
% Train it
% 

if normalizeData
  fprintf('%s: normalizing features to have mean 0 and stdev 1\n',this);
  examples = normalize(examples);
end

%labelsN = OneOfNencoding(labels); % switch labels to 1-of-N encoding
%labelsN = 2*labelsN - 1;          % now make them -1/1
%labelsN = labelsN(:,1);

if useCV
  % run over several possible values for the C parameter,
  % collecting CV error, pick a value of C and put it in the argument string

  % SVMlight default
  start = 1 / (mean(sqrt(sum(examples .* examples,2))).^2);
  
  % idea: quickly explore up and down, and then proceed on the direction where
  % the max is
  fprintf('%s: running cross-validation to determine C\n',this);
  
  current = start;
  ks = sprintf('%s -x 1 -v 0 -c %s',kernelString,num2str(current));
  modelStart = svmlearn(examples,labels,ks);
  errorStart = modelStart.loo_error;

  incrementValue = 2;

  fprintf('\tstart at C = %1.4f\terror = %1.2f\n',start,errorStart);
  
  above = current*incrementValue; below = current/incrementValue;
  tmp = 1;
  computeAbove = 1; computeBelow = 1;
  nIterations = 8;
  
  fprintf('\tabove\t\t\t\tbelow\t\t\n');
  
  while tmp <= nIterations

    if computeBelow
      ks = sprintf('%s -x 1 -v 0 -c %s',kernelString,num2str(below));
      modelBelow = svmlearn(examples,labels,ks);
      errorBelow = modelBelow.loo_error;
    end

    if computeAbove
      ks = sprintf('%s -x 1 -v 0 -c %s',kernelString,num2str(above));
      modelAbove = svmlearn(examples,labels,ks);
      errorAbove = modelAbove.loo_error;
    end
      
    fprintf('\tC = %1.6f\terror = %1.2f\tCb = %1.6f\terror = %1.2f\n',above,errorAbove,below,errorBelow);
    
    % tmp = (errorBelow == errorStart) & (errorAbove == errorStart);
    % if tmp;  above = above * incrementValue; below = below / incrementValue; end
  
    tmp = tmp + 1;;
    if (errorBelow == errorStart) & (errorAbove == errorStart)
      above = above * incrementValue;
      below = below / incrementValue;
    else
      if (errorBelow == errorStart)
	computeAbove = 0;
	below = below / incrementValue;
      elseif (errorAbove == errorStart) 
	computeBelow = 0;
	above = above * incrementValue;
      else
	% both different, stop
	break;
      end
    end
  end
  
  %% at this point, we finished because:
  %% - errorAbove and/or errorBelow are different from errorStart
  %% - we ran out of iterations
  
  if (errorBelow >= errorStart & errorAbove >= errorStart)
    % just use the default
    optimC = start;
    current      = start;
    errorCurrent = errorStart;
    moving       = 'done';
  
    fprintf('done C = %1.4f error = %1.2f\n',current,errorCurrent);  
  else
    % at least one of above and below is smaller than start
    
    if errorBelow < errorStart
      % push in the smaller direction
      changeFactor = 1/incrementValue;
      errorCurrent = errorBelow;
      current      = below;
      moving       = 'down';
      %	fprintf('\tsearch below C = %1.4f error = %1.2f\n',current,errorCurrent);
    else
      changeFactor = incrementValue;
      errorCurrent = errorAbove;
      current      = above;
      moving       = 'up';
      %	fprintf('\tsearch above C = %1.4f error = %1.2f\n',current,errorCurrent);
    end
  
    % now move in the direction found, stopping if there's no change
    fprintf('search C = %1.4f error = %1.2f moving %s\n',current,errorCurrent,moving);  
    tmp = 1;
    while tmp 
      
      next = current * changeFactor;
      
      ks = sprintf('%s -x 1 -v 0 -c %s',kernelString,num2str(next));
      modelNext = svmlearn(examples,labels,ks);
      errorNext = modelBelow.loo_error;
      
      fprintf('search C = %1.4f error = %1.2f\n',next,errorNext);
      
      if errorNext < errorCurrent
	current      = next;
	errorCurrent = errorNext;
      else
	tmp = 0;
	last      = next;
	errorLast = errorNext;
      end
    end

    % finally, do a search at regular intervals and pick the
    % smallest C that gives the minimum score
    
    nIntermediate = 5;
    if current < last
      valuesToTry   = linspace(current,last,nIntermediate);
    else
      valuesToTry   = linspace(last,current,nIntermediate);
    end
      
    fprintf('%s: searching finely between %1.4f and %1.4f\n',this,current,last);
    fprintf('\tC\tloo error\n');
    for i = 1:nIntermediate
      kernelStringIntermediate = sprintf('%s -x 1 -v 0 -c %s',kernelString,num2str(valuesToTry(i)));
      modelIntermediate{i} = svmlearn(examples,labels,kernelStringIntermediate);
      errorIntermediate(i) = modelIntermediate{i}.loo_error;    
      valueIntermediate(i) = valuesToTry(i);
      
      fprintf('\t%s\t%s\n',num2str(valueIntermediate(i)),num2str(errorIntermediate(i)));
    end
  
    [minError,minPos] = min(errorIntermediate);
    optimC = valueIntermediate(minPos);
    fprintf('\n\tusing %s %1.2f\n',num2str(optimC),errorIntermediate);
  end

  % add to kernel string
  kernelString = sprintf('%s -v 1 -c %s',kernelString,num2str(optimC));
end

%% train!
kernelString
model = svmlearn(examples,labels,kernelString);

%% output additional information

switch kernel

 case {'linear'}
  % find the hyperplane weights
  
  % training set predictions
  [ err, predictions ] = svmclassify(examples,labels,model);
    
  if 1
    % regress on the predictions to pull out weights
    if 0
      b = regress(predictions,[ones(size(examples,1),1),examples]);
      w = b(2:end); w0 = b(1);    
    else
      whos
      b = compute_multivariateRegression(predictions,examples);
      w = b(1:end-1); w0 = b(end);
    end
      
    if 0
      % check
      wpredictions = examples * w + w0;  
      [predictions,wpredictions]
      pause
    end
    
    model.extraInformation{1} = w;
    model.extraInformation{2} = w0;
  end
    
 otherwise
  % no addditional information fields
end
  
  
%
% Ancillary code 
%

%% Transform a list of labels into a "1 of N" encoding
%% (a binary matrix with as many rows as examples, and as many
%% columns as labels. The row for an example is 1 in the column
%% corresponding to its label, and 0 everywhere else. The columns
%% are in the order of the labels sorted by value.

function [labels1ofN] = OneOfNencoding(labels)

classes   = unique(labels); nClasses = length(classes);
nExamples = length(labels);

labels1ofN = zeros(nExamples,nClasses);
for c = 1:nClasses
  label           = classes(c);
  labels1ofN(:,c) = (labels == label);
end


%% Normalize each feature to have mean 0 and standard deviation 1

function [Y] = normalize(X)

[nExamples,nFeatures] = size(X);
meanX = mean(X,1);
stdvX = std(X,0,1);

Y = X -  repmat(meanX,[nExamples,1]);
Y = Y ./ repmat(stdvX,[nExamples,1]);

