function [gamma, xi, loglik] = forwards_backwards_ktattn(prior, Amat, obslikB, obslikC, obslikD)
% FORWARDS_BACKWARDS Compute the posterior probs. in an HMM using the forwards backwards algo.
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

T = size(obslikB, 2);   % T -- the length of the sequence
N = length(prior);      % N -- the number of states

scale = ones(1,T);
loglik = 0;
alpha = zeros(N,T);
gamma = zeros(N,T);
xi = zeros(N,N,T-1);

% forward
t = 1;
alpha(:,1) = prior(:) .* obslikB(:,t) .* obslikD(:,t);    % eqn (19)
[alpha(:,t), scale(t)] = normalise(alpha(:,t));
Amat2 = Amat';
for t=2:T                 % eqn (20) [note: in equations, uses t, t+1. here, uses t-1, t]
    [alpha(:,t),scale(t)] = normalise((Amat2 * (alpha(:,t-1) .* obslikC(:,t-1))) .* obslikB(:,t) .* obslikD(:,t));
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

% backward
beta = zeros(N,T); 
gamma = zeros(N,T);
beta(:,T) = ones(N,1);
gamma(:,T) = normalise(alpha(:,T) .* beta(:,T));
t=T;
for t=T-1:-1:1
    b = beta(:,t+1) .* obslikB(:,t+1) .* obslikC(:,t) .* obslikD(:,t+1);   % eqn (25)
    beta(:,t) = normalise((Amat * b));  % eqn (25)
    gamma(:,t) = normalise(alpha(:,t) .* beta(:,t));  % eqn (27)
    xi(:,:,t) = normalise((Amat .* (alpha(:,t) * b')));  % eqn (37)
end

