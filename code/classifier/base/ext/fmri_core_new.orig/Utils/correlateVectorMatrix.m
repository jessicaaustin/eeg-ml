% Computes the correlation of a query vector of observations
% with the columns of a matrix.
%
% In:
% - vector (column vector, rows are observations)
% - matrix (where each column is a variable, and rows are observations)
%
% Out:
% - a row vector with correlation between the probe vector and each
% row in the matrix
%
% Dependencies:

% History
% - 30 May 03 - fp - created
%

function [result] = correlateVectorMatrix( varargin )
  
  l = length(varargin);
  if l < 1
    fprintf(1,'syntax: correlateVectorMatrix(v,m)\n');
    return;
  end
  v = varargin{1};
  m = varargin{2};
  clear varargin;
  
  [h,w] = size(v); % make sure it is a column vector
  if h == 1; v = v'; end
  
  %% Figure out things
  
  [nEx,nVars] = size(m);
  
%  vm         = v * ones(1,nVars);
  vm         = repmat(v,1,nVars);
  meanMatrix = mean(m,1);
  meanVector = mean(vm,1);
  stdMatrix  = std(m,1,1);
  stdVector  = std(vm,1,1);
  tmp        = m .* vm;
  result     = (mean(tmp,1) - (meanMatrix .* meanVector)) ./ (stdMatrix .* stdVector);
  
 
function [] = testThis()

a = randn(100,3);
m = zeros(100,9);

m(:,1) = a(:,1);
m(:,2) = a(:,1) + randn(100,1)/0.01;
m(:,3) = -1*a(:,1);
m(:,4) = a(:,2);
m(:,5) = a(:,2) + randn(100,1)/0.01;
m(:,6) = -1*a(:,2);
m(:,7) = a(:,3);
m(:,8) = a(:,3) + randn(100,1)/0.01;
m(:,9) = -1*a(:,3);

