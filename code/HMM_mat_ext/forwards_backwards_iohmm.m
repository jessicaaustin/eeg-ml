function [gamma, xi, loglik] = forwards_backwards_iohmm(prior, transmat, obslik, maximize)
% FORWARDS_BACKWARDS_IOHMM Compute the posterior probs. in an HMM using the forwards backwards algo.
% [gamma, xi, loglik] = forwards_backwards(prior, transmat, obslik, maximize)
% Use obslik = mk_dhmm_obs_lik(data, b) or obslik = mk_ghmm_obs_lik(data, mu, sigma) first.
%
% Inputs:
% PRIOR(I) = Pr(Q(1) = I)
% TRANSMAT(I,J) = Pr(Q(T+1)=J | Q(T)=I)
% OBSLIK(I,T) = Pr(Y(T) | Q(T)=I)
% maximize is optional; if 1, we do max-product (as in Viterbi) instead of sum-product
%
% Outputs:
% gamma(i,t) = Pr(X(t)=i | O(1:T))
% xi(i,j,t)  = Pr(X(t)=i, X(t+1)=j | O(1:T)) t <= T-1

if nargin < 4, maximize = 0; end

T = size(obslik, 2);
Q = length(prior);

scale = ones(1,T);
loglik = 0; 
alpha = zeros(Q,T); 
gamma = zeros(Q,T);
xi = zeros(Q,Q,T-1);

% FORWARD (aka Filtering) -- find p(ht|v1:t)   [Barber §23.2.2]
% α(ht) = p(ht, v1:t)

t = 1;  % α(h1) = p(h1, v1) = p(v1|h1)p(h1)     [Barber 23.2.10]
alpha(:,1) = prior(:) .* obslik(:,t);
[alpha(:,t), scale(t)] = normalise(alpha(:,t));
transmat2 = transmat';
for t=2:T
  if maximize
    A = repmat(alpha(:,t-1), [1 Q]);
    m = max(transmat .* A, [], 1);
    [alpha(:,t),scale(t)] = normalise(m(:) .* obslik(:,t));
  else  %                              p(ht|ht−1)    α(ht−1)          p(vt|ht)    [Barber 23.2.9]
    [alpha(:,t),scale(t)] = normalise((transmat2 * alpha(:,t-1)) .* obslik(:,t));
  end
  if (scale(t) == 0) | isnan(scale(t)) | ~isreal(scale(t)) 
    fprintf('scale(%d)=%5.3f\n', t, scale(t))
    keyboard
  end
end
if any(scale==0)
  loglik = -inf;
else
  loglik = sum(log(scale));
end

% BACKWARD (aka Parallel smoothing) -- find p(ht|v1:T)  [Barber 23.2.3]
% β(ht) ≡ p(vt+1:T |ht)

beta = zeros(Q,T); % beta(i,t)  = Pr(O(t+1:T) | X(t)=i)
gamma = zeros(Q,T);
beta(:,T) = ones(Q,1);
gamma(:,T) = normalise(alpha(:,T) .* beta(:,T));
t=T;
for t=T-1:-1:1
  %      β(ht+1)     p(vt+1|ht+1)     [Barber 23.2.18]
  b = beta(:,t+1) .* obslik(:,t+1); 
  if maximize
    B = repmat(b(:)', Q, 1);
    beta(:,t) = normalise(max(transmat .* B, [], 2));
  else               %    p(ht+1|ht)                    [Barber 23.2.18]
    beta(:,t) = normalise((transmat * b));
  end
  gamma(:,t) = normalise(alpha(:,t) .* beta(:,t));   % [Barber 23.2.19]
  xi(:,:,t) = normalise((transmat .* (alpha(:,t) * b')));
end


