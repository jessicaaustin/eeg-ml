% Train a K class logistic regression classifier. Differs from
% other versions in the core in that it allows specification
% of the optimization algorithm to use and regularization method.
%
% In:
% - examples
% - labels
% - optional parameters (add arguments (examples,labels,<parameterName>,<parameterValue>,...)
%   - 'regularization' - which regularization penalty to use
%                        values: 'none'|'L1'|'L2' (default)
%   - 'lambda' - the value of the tradeoff parameter for regularization (default 1)
%
%   - 'optimization' - two parameters:
%     - which optimization approach to use (default: 'gradientDescent')
%     - a cell array of parameters for the optimization algorithm
%
% Out:
% - a models structure containing the classifier
%   (please see comments in trainClassifier for more details)
%
% Dependencies:
%
% History:
% - 2007 Jun 6 - fpereira@cs.cmu.edu
%
% Examples:

% [models] = classifierLogisticRegressionMulti(examples,labels);
% [models] = classifierLogisticRegressionMulti(examples,labels,'regularization',{'L2',1});
% [models] = classifierLogisticRegressionMulti(examples,labels,'regularization',{'L1',1});
%
% Notes:
%
% 


function [models] = classifierLogisticRegressionMulti( varargin )

%
% process arguments
%

this = 'classifierLogisticRegression3';

if nargin < 2; eval(sprintf('help %s;',this)); return; else

  examples  = varargin{1};
  labels    = varargin{2};

  % defaults for options
  regularization = 'L2';
  lambda = 1;
  regularizationParameters = {regularization,lambda};
  optimization = 'gradientDescent';
  step = 0.01;
  optimizationParameters = {optimization,step}; % step size
  voxelsToNeighbours = []; numberOfNeighbours = [];
  
  if nargin > 2
    % there are additional arguments to process
  
    idx = 3;
    while idx <= nargin
      argName = varargin{idx}; idx = idx + 1;

      switch argName
       case {'regularization'}
	regularizationParameters = varargin{idx}; idx = idx + 1;
	regularization = regularizationParameters{1};
	switch regularization
	 case {'L1','L2'}
	  lambda = regularizationParameters{2};
	end
       
       case {'optimization'}
	optimizationParameters = varargin{idx}; idx = idx + 1;
	optimization = optimizationParameters{1};
	switch optimization
	 case {'gradientDescent'}
	  step = optimizationParameters{2};
	end
	
       case {'neighbourhood'}
	voxelsToNeighbours = varargin{idx}; idx = idx + 1;
	numberOfNeighbours = varargin{idx}; idx = idx + 1;
       otherwise
	fprintf('%s: unknown argument %s\n',this,argName); pause
      end 
    end
  end

end

%% compute information

sortedLabelValues     = sort(unique(labels)); conditions = sortedLabelValues;
nClasses              = length(sortedLabelValues);
[nExamples,nFeatures] = size(examples);


%% find examples in each class and create a binary indicator mask

delta = zeros(nExamples,nClasses);

for c = 1:nClasses
  label  = sortedLabelValues(c);  
  indices{c}  = find(labels == label);
  nperclass(c) = length(indices{c});

  delta(indices{c},c) = 1;
end


%% compute class mean examples

meanExamples = zeros(nClasses,nFeatures);
meanLabels   = zeros(nClasses,nClasses);

for c = 1:nClasses
  meanExamples(c,:) = mean(examples(indices{c},:),1);
  meanLabels(c,1)   = sortedLabelValues(c);
end  


%
% Train classifier
%

X = [ones(nExamples,1),examples];
nFeatures = nFeatures + 1;
Y = labels;
%W = randn(nFeatures,nClasses);
W = zeros(nFeatures,nClasses);

%W = repmat(0.1,nFeatures,nConds);


switch optimization

 case {'gradientDescent'}
  % all regularizations share the same loop

  fprintf('%s: with %s regularization, alpha=%d\n',this,regularization,lambda);
  
  switch regularization
   case {'L1'}
    [W,obj] = internal_gradientdescent('internal_error_logregL1',{step},{W,X,Y,lambda,delta,conditions});
   case {'L2'}
    [W,obj] = internal_gradientdescent('internal_error_logregL2',{step},{W,X,Y,lambda,delta,conditions}); 
  end
end 


%% Pack the results into a cell array

% Discriminative model - a cell array of sets of weights
models{nClasses+1} = W;

% Generative model - nothing here

% Training Set information
trainingSetInfo.nExamples         = nExamples;
trainingSetInfo.nFeatures         = nFeatures;
trainingSetInfo.nLabels           = nClasses;
trainingSetInfo.nClasses          = nClasses; % same thing as labels
trainingSetInfo.sortedLabelValues = sortedLabelValues;
trainingSetInfo.classPriors       = nperclass / sum(nperclass);
models{nClasses+2} = trainingSetInfo;

% Extra information (depends on classifier)

models{nClasses+3}{1} = regularization;
models{nClasses+3}{2} = regularizationParameters;
models{nClasses+3}{3} = optimization;
models{nClasses+3}{4} = optimizationParameters;
models{nClasses+3}{5} = voxelsToNeighbours;
models{nClasses+3}{6} = numberOfNeighbours;
models{nClasses+3}{7} = meanExamples;

%
% Ancillary functions
%

function [W,f] = internal_gradientdescent(varargin)

functionToCall     = varargin{1};
descentParameters  = varargin{2};
functionParameters = varargin{3};

% initial settings
W = functionParameters{1}; 
[nFeatures,nConds] = size(W);
W(:,nConds) = repmat(0,nFeatures,1); % make sure

f = -Inf; df = repmat(-Inf,size(W));

% create a function call that can be eval()-ed
argString = 'W';
for n = 2:length(functionParameters);
  argString = [argString,',',sprintf('functionParameters{%d}',n)];
end
callString = sprintf('[f,df] = %s(%s);',functionToCall,argString);

f_previous  = f;
W_previous  = W;
df_previous = df;

%% start iterating

step      = 0.01;
threshold = 0.00001;
if ~isempty(descentParameters)
  ldp = length(descentParameters);
  step = descentParameters{1};
  if ldp > 1; threshold = descentParameters{2}; end
end  
  

eval(callString); % sets f and df

%fprintf('0: %s\n',num2str(f));

iteration = 1;
warmedup  = 0;

while 1
  
  W = W + step * df;
  W(:,nConds) = repmat(0,nFeatures,1);
  
  % call error function to give us new:
  % - f  - function value at x
  % - df - gradient at x
  eval(callString);
  
  if ~rem(iteration,100); fprintf('%d: %1.8f\n',iteration,f); end
  
  if f < f_previous
    % went past, reset to previous values and retry lower step
    fprintf('gradientdescent: %d: objective decrease, decimating step and backtracking\n',iteration);
    step = step / 10;

    W  = W_previous;
    f  = f_previous;
    df = df_previous;
  
  else
    change = (f-f_previous)/abs(f);
    
    if warmedup
      
      if isnan(f) | (abs(f)<(eps*10^4)) | (change < threshold)
	% no tangible progress, jump out
	if isnan(f)
	  fprintf('gradientdescent: f is NaN, stopping\n');
	else
	  fprintf('gradientdescent: %d: change is %s, stopping\n',iteration,num2str(change));
	end
	break;
      else
	% keep going and accelerate
	step = 1.1 * step;
      end
    
    else
      warmedup = 1;
    end
    
  end

  W_previous  = W;
  f_previous  = f;
  df_previous = df; 
  
  iteration = iteration + 1;
end
fprintf('gradientdescent: %d: converged at %1.8f\n',iteration,f);


% Computer logistic regression objective (L2 penalty) and gradient
%
% In:
% - W - weights - #features (includes a column of ones) x #conditions
% - X - data - #examples x #features (first column is ones)
% - Y - labels - #examples x 1
% - lambda - tradeoff parameter between likelihood and L2 penalty
% - conditions - same as unique(Y), there to save computational effort
%
% Out:
% - f  - function value at random point <w>
% - df - gradient at random point <w>

% notes:
% - X already includes a column of 1s
% - W(1,:) is the bias term for each class
% - L2 penalty is on W(2:end,:)
% - this works *as long as the last column of W is 0*

function [f,df] = internal_error_logregL2(W,X,Y,lambda,delta,conditions)

[nExamples,nFeatures] = size(X);
nConds = length(conditions);

df = zeros(nFeatures,nConds);

%% compute P(Y=k|X) for all X (rows) and K (cols) at once  
  
% exp(X * weights), for all K conditions (columns) over all examples (rows)  
tmp  = exp(X*W);
dsum = sum(tmp,2); % the denominator sum over all exp(X*weights)

% compute P(Y=k|X) for each k, over all examples at once,
% works *as long as the last column of W is 0*
py_k = tmp ./ repmat(dsum,1,nConds); clear tmp;


%% compute gradient

% compute all the weight updates for condition <k>
for k = 1:nConds - 1    
  errork = delta(:,k) - py_k(:,k);
  df(:,k) = sum( X .* repmat(errork,1,nFeatures),1)';
end

df = df - lambda*W;
  
%% compute log likelihood and function value

L_examples = log( sum(delta .* py_k,2) ); % all 0 except right class
L          = sum(L_examples,1); % over all data points

%f = L - 0.5*lambda*norm(W,'fro')^2;
f = L - 0.5*lambda*norm(W(2:end,1:(end-1)),'fro')^2;


% Computer logistic regression objective (L1 penalty) and gradient
%
% In:
% - W - weights - #features (includes a column of ones) x #conditions
% - X - data - #examples x #features (first column is ones)
% - Y - labels - #examples x 1
% - lambda - tradeoff parameter between likelihood and L2 penalty
% - conditions - same as unique(Y), there to save computational effort
%
% Out:
% - f  - function value at random point <w>
% - df - gradient at random point <w>

% notes:
% - a - is the parameter that regulates how good/slow the
% approximation to L1 norm and its derivative is
% - X already includes a column of 1s
% - W(1,:) is the bias term for each class
% - L2 penalty is on W(2:end,:)
% - this works *as long as the last column of W is 0*

function [f,df] = internal_error_logregL1(W,X,Y,lambda,delta,conditions)

a = 4; % higher a -> better, slower approximation to L1

[nExamples,nFeatures] = size(X);
nConds = length(conditions);

df = zeros(nFeatures,nConds);

%% compute P(Y=k|X) for all X (rows) and K (cols) at once  
  
% exp(X * weights), for all K conditions (columns) over all examples (rows)  
tmp  = exp(X*W);
dsum = sum(tmp,2); % the denominator sum over all exp(X*weights)

% compute P(Y=k|X) for each k, over all examples at once,
% works *as long as the last column of W is 0*
py_k = tmp ./ repmat(dsum,1,nConds); clear tmp;


%% compute gradient

% compute all the weight updates for condition <k>
for k = 1:nConds - 1    
  errork = delta(:,k) - py_k(:,k);
  df(:,k) = sum( X .* repmat(errork,1,nFeatures),1)';
end

df = df - lambda*((1-exp(-2*a*W))./(1+exp(-2*a*W)));

%% compute log likelihood and function value

L_examples = log( sum(delta .* py_k,2) ); % all 0 except right class
L          = sum(L_examples,1); % over all data points

%f = L - 0.5*lambda*norm(W,'fro')^2;
f = L - 0.5*lambda*(1/a)*sum(sum(log(exp(a*W)+exp(-a*W))));


