%
% Standardizes each example
%
% In:
% - examples
% - labels
% - method to use
%   - 'mean0' - the mean of each example becomes 0
%   - 'std1'  - the stdv of each example becomes 1
%   - 'mean0std1' - (default) - mean0 + std1
%   - 'correlation' - same thing, but also multiply each example by the normalized mean class example
%   - 'centerFeatures' - the mean of each feature becomes 0
%
% Out:
% - normalized examples
% - labels
% - the class means of the normalized examples
%
% Examples:
% [examples,labels,classMeanExamples] = normalizeELE_standardizeExamples(examples,labels,'mean0std1');
% (makes each example have mean 0 and standard deviation 1)
%
% History:
% - 2005 Aug 12 - fpereira - created
%

function [examples,labels,classMeanExamples] = normalizeELE_standardizeExamples(varargin)

%
% Process parameters
%

this = 'normalizeELE_standardizeExamples';
l = length(varargin);
if l < 3
  eval(sprintf('help this;'));return;
else
  examples = varargin{1};
  labels   = varargin{2};
  method = varargin{3};
end

[nExamples,nFeatures] = size(examples);

conditions = unique(labels); nConds = length(conditions);
nLabels = size(labels,2); nFeatures = size(examples,2);

classMeanExamples = zeros(nConds,nFeatures);

%% find images of each class, also mean and stdv class image

for c = 1:nConds
  cond  = conditions(c);  
  examplesCond{c}  = find(labels == cond);
  nExamplesCond(c) = length(examplesCond{c});

  classMeanExamples(c,:) = mean(examples(examplesCond{c},:),1);
end  
featureMeans = mean(examples,1);
featureStdvs = std(examples,0,1);

%
% normalize examples, if required
%

switch method

 case {'centerFeatures'}
  examples = examples - repmat(featureMeans,nExamples,1);
 
 otherwise
  % case {'mean0','std1','mean0std1','correlation','mean0std1_forCorrelation','norm1'}
  exampleMeans = mean(examples,2);
  exampleStdvs = std(examples,0,2);
  
  switch method     
   case {'mean0'}
    examples = examples -  repmat(exampleMeans,1,nFeatures);         
   case {'std1'}
    examples = examples ./ repmat(exampleStdvs,1,nFeatures);
   case {'mean0std1'}
    examples = examples -  repmat(exampleMeans,1,nFeatures);
    examples = examples ./ repmat(exampleStdvs,1,nFeatures);
   case {'mean0std1_forCorrelation'}
    % MATLAB correlation function divides by std estimate with
    % denominator N rather than (N-1)
    exampleStdvs = std(examples,1,2);
    examples = examples -  repmat(exampleMeans,1,nFeatures);
    examples = examples ./ repmat(exampleStdvs,1,nFeatures);
   case {'norm1'}
    exampleNorms = sqrt(sum(examples .* examples,2));
    examples = examples ./ repmat(exampleNorms,1,nFeatures);
    
   case {'correlation'}
    examples = examples -  repmat(exampleMeans,1,nFeatures);
    examples = examples ./ repmat(exampleStdvs,1,nFeatures);
    
    % mean 0 std1 class mean examples (needed for correlation)
    cmeMeans = mean(classMeanExamples,2);
    cmeStdvs = std(classMeanExamples,0,2);
    classMeanExamples = classMeanExamples -  repmat(cmeMeans,1,nFeatures);
    classMeanExamples = classMeanExamples ./ repmat(cmeStdvs,1,nFeatures);
  
  end
end

%
% Compute
%

switch method

  
 case {'mean0','std1','mean0std1','centerFeatures','mean0std1_forCorrelation'}
  % all done
  
 case {'correlation'}
  
  %
  % multiply each example by its class mean, component-wise
  %
  
  for c = 1:nConds
    cond  = conditions(c);  
    tmp   = repmat(classMeanExamples(c,:),nExamplesCond(c),1);
    examples(examplesCond{c},:) = examples(examplesCond{c},:) .* tmp;
  end
  
 otherwise
end
