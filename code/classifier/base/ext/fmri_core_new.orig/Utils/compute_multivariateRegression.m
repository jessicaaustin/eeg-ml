% Computes a multivariate multiple linear regression, that is
%
% Input:
% - a matrix Y where each column is the *dependent* variable on a
%   multiple linear regression
% - a matrix X where each column is an *independent* variable to be
%   used on the multiple linear regression
%
% in a nutshell, each column of Y is predicted as a linear
% combination of all the columns of X.
%
% Out:
% - betas
%     a |X|+1 x |Y| matrix. Each column contains the regression coefficients
%     of the corresponding column of Y on all the X (first |X| rows) and 1 (last row).
% - residuals
%     = Y - X Betas
%
% History:
% 30 Dec 2004 - fpereira - created
                                                                                                  
function [betas,residuals] = compute_multivariateRegression( varargin )
                                                                                                  
l = length(varargin);
  
if l < 1
  help compute_multivariateRegression; return;
else
  Y = varargin{1};
  [nExamples,nVars] = size(Y);
    
  % take the X provided or use 1:N as the X
  if l > 1;
    X = [varargin{2},ones(nExamples,1)];
  else
    X = [(1:nExamples)',ones(nExamples,1)];
  end
end

% used to figure out whether the usual method, which requires a
% matrix inversion, is feasible, as rcond(matrix) tells us that.
tmp = X'*X;
r1  = rank(tmp); r = 1;
r2  = rcond(tmp);

%% Compute actual 

if rcond(tmp) > 0.00001

  fprintf('compute_multivariateRegression: using regular algorithm\n');

  % we can invert it, so use the normal formula
  betas = inv(X'*X) * X' * Y;

else

  fprintf('compute_multivariateRegression: rank(X''X)=%d, using SVD(X)\n',r);
  
  % Y = X * L is the regression problem, we want L, but X'X is not invertible
  %
  % X = ABC' (SVD) so
  % Y = ABC' L
  % Y = AB(C'L) = AB (L'C)'
  % which is a regression Y on AB with coefficients M = (L'C)',
  % doable because AB can be inverted. Given that, M can be determined,
  % and M = (L'C)' solved for L, yielding L = CM
  %
  % Computing the regression on AB is also simplified, as
  % inv( (AB)'AB ) = inv( B'A'A B ) = inv( B'B ), as A'A = I
  % inv(B'B) (AB)' = B^-1 B'^-1 B' A' = B^-1 A'
  
  [A,B,C] = compute_fastSVD(X);
  
  % reduce this
  r = min(size(A),size(C));
  d = diag(B);
  d = d / sum(d);
  d = cumsum(d);
  tmp = find(d >= 0.99); pos = tmp(1);
  A = A(:,1:pos);
  C = C(:,1:pos);
  B = B(1:pos,1:pos);

  % compute 
  M = inv(B')*A'  * Y;
  betas = C*M;
end
 
residuals = Y - X * betas;











function [] = testThis()

% simple case

nVars = 10;
nTime = 100;

X = rand(nTime,nVars)*2;
r = randperm(nVars);
Y = X * r' + randn(nTime,1);

% similar, but X now has 2 copies of every column
X = [X,X];

% similar, but the 2 copies have independent noise

X  = rand(nTime,nVars)*2;
X1 = X + randn(nTime,nVars);
X2 = X + randn(nTime,nVars);

r = randperm(nVars);
Y = X1*r' + X2*r' + randn(nTime,1);


% same thing but several Y

nY = 10;
Y  = zeros(nTime,nY);
r  = zeros(nY,nVars);

for y = 1:nY  
  r(y,:) = randperm(nVars);
  Y(:,y) = X1 * r(y,:)' + X2 * r(y,:)' + randn(nTime,1);
end





% fewer time points than variables

nVars = 500;
nTime = 100;

X  = rand(nTime,nVars)*2;
r1 = randperm(nVars);

Y = X * r1 + randn(nTime,1);

% same thing but several Y

nY = 10;
Y  = zeros(nTime,nY);
r  = zeros(nY,nVars);

for y = 1:nY  
  r(y,:) = randperm(nVars);
  Y(:,y) = X * r(y,:)' + randn(nTime,1);
end



%Y(:,1) = 1*X + 4 + randn(nTime,1);
%Y(:,2) = 2*X + 5 + randn(nTime,1);
%Y(:,3) = 3*X + 6 + randn(nTime,1);

