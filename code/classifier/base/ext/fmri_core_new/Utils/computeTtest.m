function [h,pvals,ci,stats] = computeTtest(X,m,alpha,tail)
%TTEST  Hypothesis test: Compares the sample average to a constant.
%
%   Modified version of the matlab ttest code that does many tests
%   simultaneously - each column of X is a sample.
%
%   [H,P,CI,STATS] = TTEST(X,M,ALPHA,TAIL) performs a T-test to determine
%   if a sample from a normal distribution (in X) could have mean M.
%   M = 0, ALPHA = 0.05 and TAIL = 0 by default.
%
%   The Null hypothesis is: "mean is equal to M".
%   For TAIL=0,  alternative: "mean is not M".
%   For TAIL=1,  alternative: "mean is greater than M"
%   For TAIL=-1, alternative: "mean is less than M"
%   TAIL = 0 by default.
%
%   ALPHA is desired significance level. 
%   P is the p-value, or the probability of observing the given result
%     by chance given that the null hypothesis is true. Small values
%     of P cast doubt on the validity of the null hypothesis.
%   CI is a confidence interval for the true mean.  Its confidence
%     level is 1-ALPHA.
%   STATS is a structure with two elements named 'tstat' (the value
%     of the test statistic) and 'df' (its degrees of freedom).
%
%   H=0 => "Do not reject null hypothesis at significance level of alpha."
%   H=1 => "Reject null hypothesis at significance level of alpha."

%   References:
%      [1] E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, 1970, page 206. 

%   Copyright 1993-2000 The MathWorks, Inc. 
%   $Revision: 2.11 $  $Date: 2000/05/26 18:53:50 $

if nargin < 1, 
    error('Requires at least one input argument.'); 
end

[n1 ntests] = size(X);
%x = x(~isnan(x));

if nargin < 2
    m = 0;
end
 
if nargin < 4, 
    tail = 0; 
end 

if nargin < 3, 
    alpha = 0.05; 
end 

if (prod(size(alpha))>1), error('ALPHA must be a scalar.'); end
if (alpha<=0 | alpha>=1), error('ALPHA must be between 0 and 1'); end

samplesize  = n1;
xmeans = mean(X);
sers = std(X,0,1) ./ sqrt(samplesize);
tvals = (xmeans - m) ./ sers;
indices = find(sers==0);
%if length(indices)
%  sers
%  pause
%end

if (nargout > 3), stats = struct('tstat', tvals, 'df', samplesize-1); end;
pvals = tcdf(tvals,samplesize - 1);

% the p-value just found is for the  tail = -1 test
crit = tinv(1 - alpha,samplesize - 1) .* sers;

if (tail == 0)
    pvals = 2 * min(pvals, 1-pvals);
    crits = tinv((1 - alpha / 2),samplesize - 1) .* sers;
    if (nargout>2), ci = [(xmeans - crits); (xmeans + crits)]; end
else
    crits = tinv(1 - alpha,samplesize - 1) .* sers;
    if tail == 1
        pvals = 1 - pvals;
        if (nargout>2), ci = [(xmeans - crits); repmat(Inf,1,ntests)]; end
    else
        if (nargout>2), ci = [repmat(-Inf,1,ntests); (xmean + crit)]; end
    end
end

% Determine if the actual significance exceeds the desired significance
h = (pvals <= alpha);



function [] = testThis()

Y = randn(100,5);
Y = Y + ones(100,1) * [0 0.25 0.5 0.75 1];

for v = 1:5
  [h1{1}(v),pv1{1}(v)] = ttest(Y(:,v),0,0.05,-1);
  [h1{2}(v),pv1{2}(v)] = ttest(Y(:,v),0,0.05,0);
  [h1{3}(v),pv1{3}(v)] = ttest(Y(:,v),0,0.05,1);
end

[h2{1},pv2{1}] = ttestparallel(Y,0,0.05,-1);
[h2{2},pv2{2}] = ttestparallel(Y,0,0.05,0);
[h2{3},pv2{3}] = ttestparallel(Y,0,0.05,1);
