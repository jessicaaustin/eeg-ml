%
% Code for computing various distances measures (some are metrics,
% some aren't) between examples in two sets.
%
% In:
% - distance measure
%     euclidean, absoluteValue, correlation, cosine
%     also correlationAsDistance (convert correlation range from [-1,1] to [0,1], where 0
%     corresponds to 1 (perfectly correlated is closest) and 1 to -1)
%
% - example set 1, a <#examples 1> x <#features> matrix
% - optional (default is to use example set 1)
%   - example set 2, a <#examples 2> x <#features> matrix
%
% Out:
% - a matrix with <#examples 1> x <#examples 2> values,
%   where entry <i,j> has the distance measure between
%   the i-th example of set 1 and the j-th example of set 2
%
% Examples:
% [measureValues] = compute_distanceMeasures('euclidean',examples1);
% [measureValues] = compute_distanceMeasures('correlation',examples1,examples2);
%
% Notes:
%
% History:
% - 2005 Aug 10 - fpereira - added correlationAsDistance
% - 2005 Jul 22 - fpereira - created
%

function [measureValues] = compute_distanceMeasures( varargin )

%
% Process parameters
%

this = 'compute_distanceMeasures';
l    = length(varargin);
if l < 2; eval(sprintf('help %s',this)); return; end

measure   = varargin{1};
examplesA = varargin{2};
if l > 2; examplesB = varargin{3}; else; examplesB = examplesA; end

[nA,nFeaturesA] = size(examplesA);
[nB,nFeaturesB] = size(examplesB);

if nFeaturesA ~= nFeaturesB
  eval(sprintf('help %s',this)); return;
else nFeatures = nFeaturesA; end

%
% Compute
%

switch measure

 case {'euclidean'}
  measureValues = zeros(nA,nB);  

  for n = 1:nB
    measureValues(:,n) = sqrt(sum((examplesA-repmat(examplesB(n,:),nA,1)).^2,2));
  end

 case {'absoluteValue'}
  measureValues = zeros(nA,nB);  
  for n = 1:nB
    measureValues(:,n) = sum(abs(examplesA-repmat(examplesB(n,:),nA,1)),2);
  end

 case {'cosine'}
  [examplesA,dummy] = normalizeELE_standardizeExamples(examplesA,ones(nA,1),'norm1');
  [examplesB,dummy] = normalizeELE_standardizeExamples(examplesB,ones(nB,1),'norm1');
  measureValues = examplesA * examplesB';
  
 case {'correlation','correlationAsDistance'}
  [examplesA,dummy] = normalizeELE_standardizeExamples(examplesA,ones(nA,1),'mean0std1_forCorrelation');
  [examplesB,dummy] = normalizeELE_standardizeExamples(examplesB,ones(nB,1),'mean0std1_forCorrelation');
  measureValues = (examplesA * examplesB')/nFeatures; clear tmp;

  if isequal(measure,'correlationAsDistance')
    % convert correlation range from [-1,1] to [0,1], where 0
    % corresponds to 1 (perfectly correlated is closest) and 1 to -1
    measureValues = -0.5*measureValues + 0.5;
  end
  
 otherwise
  fprintf('%s: measure %s is not supported\n',this,measure); pause; return
    
end
