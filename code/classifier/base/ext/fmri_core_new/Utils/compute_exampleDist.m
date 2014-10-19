% Computes all pairwise distances between two sets of examples
%
% Right now, the code loops over the first set of examples
% computing distances between each of those and the second set.
%
% Properties:
% - no matter what distance is used, the result is such that
% the smaller the distance between two vectors, the closer they
% are
% - distance = 0 means the two vectors are the same
% - distance >= 0 otherwise
%
% Exceptions:
% - realCorrelation is correlation, hence lies in [-1,1]
% (the 'correlation' distance measure maps [-1,1] to the
% [0,1] interval, with 0 corresponding to 1 (perfect correlation)
% and 1 corresponding to -1 (perfect anticorrelation)

% In:
% - two matrices with dimensions #examples * #dimensions
%   (i.e. each example is a row)
% - optionally, one or more parameters to specify a distance measure
%
% Out:
% - Distances - #examples1 * #examples2 where Distances(i,j) is
%               the requested distance between example <i> in the
%               first matrix and example <j> in the second
%
%
% History:
% - 11 Mar 05 - fpereira - adapted from the earlier version
% - 27 Mar 03 - fpereira - created, following discussion with and
% implementation by pbennett
%

function [Distances,DimDistances] = compute_exampleDist( varargin )

%
% Process parameters
%

l = length(varargin);
  
if l < 2; help compute_exampleDist; return; else;
  
  E1 = varargin{1}; varargin{1} = [];
  E2 = varargin{2}; varargin{2} = [];
  [numExamples1,numFeatures1] = size(E1);
  [numExamples2,numFeatures2] = size(E2);
  
  measure  = 'euclidean'; % by default
  weights = ones(numFeatures1,1);
  if l > 2
    measure = varargin{3};
    switch measure
     case {'euclidean','statistical','hamming','cosine'}
      if l > 3; weights = varargin{4}; end
     case {'correlation','realCorrelation'}
      if l > 3; lag = varargin{4}; else; lag = 0; end
     otherwise
      fprintf('compute_exampleDist: error: measure %d is not supported\n');pause;return;
    end
  end
end

% a few checks and calculations

if numFeatures1 ~= numFeatures2
  fprintf('exampleDist: error: examples should have the same number of features in both sets\n');pause;return;
end
numFeatures = numFeatures1;

if length(weights) ~= numFeatures
  fprintf('exampleDist: error: there should be as many weights as dimensions\n');pause;return;
end

%
% Do actual distance computation
%

% Distances between points in rows and points in cols
Distances    = zeros(numExamples1,numExamples2);

switch measure
  
  
 case {'euclidean'}  
  for e=1:numExamples1
    
    E           = repmat(E1(e,:),numExamples2,1);
    Differences = E - E2; % each row is the difference between example e and an example in E2
    clear('E');
    tmp = sum(Differences .^ 2,2); clear('Differences');
			  
    Distances(e,:) = sqrt(tmp)';
    %       Differences = Differences.^2; % slow (Mat Mult)?
    %       Distances(e,:)      = (sqrt(sum((Differences .* Differences,2)))';

  end

  
 case {'correlation','realCorrelation'}

  % samples are voxels
  % variables are the different images    
  if isequal(E1,E2)
    Distances = corrcoef(E1'); % #examples * #examples
    DimDistances = [];
  else            
    for e1=1:numExamples1
      % Compute the correlation between the #features-dimensional
      % vector e1 and each row vector in E2	
      %
      % results is a row vector with #observations dimensions
      results = correlateVectorMatrix(E1(e1,:)',E2');
      Distances(e1,:) = results;
    end
  end
  
  if isequal(measure,'correlation')
    % finally invert the scale so that a smaller value means
    % that examples are closer. Final scale is in the [0,1] range.
    % with 0 being perfectly correlated, 0.5 uncorrelated, 1 anticorrelated
    Distances = -0.5 * Distances + 0.5;
  else
    % do nothing
  end

  
 case {'cosine'} 
  
  for e=1:numExamples1
    
    testNorms = sqrt(sum( E2 .* E2 ,2));
    trainNorm = norm(E1(e,:)) * testNorms;
    
    dotProds = E2 * E1(e,:)';
    dotProds = dotProds ./ trainNorm;
    
    % larger is closer, so invert scale
    Distances(e,:) = -1 * dotProds';
  end

  
 case {'hamming'}  

  for e=1:numExamples1
    
    E           = repmat(E1(e,:),numExamples2,1);
    Differences = E - E2; % each row is difference between example e and all examples in E2
    clear('E');
    Distances(e,:) = sum(abs(Differences * WeightMatrix),2)';
    
  end
  
  
 case {'statistical'}
  % same as euclidean, but each dimension is standardized by
  % dividing it by the standard deviation over the training set
  stds = std(E2,0,1);
  %    stds = ones(1,size(E2,2));
  tmp  = repmat(stds,numExamples2,1);
  E2   = E2 .* tmp;
  
  for e=1:numExamples1
    
    E           = repmat(E1(e,:)./stds, numExamples2, 1);
    Differences = E - E2; % each row is the difference between example e and an example in E2
    clear('E');
    tmp = sum(Differences .^ 2,2); clear('Differences');
    
    Distances(e,:) = sqrt(tmp)';
    %       Differences = Differences.^2; % slow (Mat Mult)?
    %       Distances(e,:)      = (sqrt(sum((Differences .* Differences,2)))';
  end
    
 otherwise
  
end
