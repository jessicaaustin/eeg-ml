function [LL, prior, Amat, Bmat, Cmat, Dmat, gamma] = learn_dhmm_ktattn(data, prior, Amat, Bmat, Cmat, Dmat, max_iter, thresh, verbose)
% LEARN_HMM Find the ML parameters of an HMM with discrete outputs using EM.
%
% [LL, PRIOR, TRANSMAT, OBSMAT] = LEARN_HMM(DATA, PRIOR0, TRANSMAT0, OBSMAT0)
% computes maximum likelihood estimates of the following parameters,
% where, for each time t, Q(t) is the hidden state, and
% Y(t) is the observation
%   prior(i) = Pr(Q(1) = i)
%   transmat(i,j) = Pr(Q(t+1)=j | Q(t)=i)
%   obsmat(i,o) = Pr(Y(t)=o | Q(t)=i)
% It uses PRIOR0 as the initial estimate of PRIOR, etc.
%
% Row l of DATA is the observation sequence for example l. If the sequences are of
% different lengths, you can pass in a cell array, so DATA{l} is a vector.
% If there is only one sequence, the estimate of prior will be poor.
% If all the sequences are of length 1, transmat cannot be estimated.
%
% LL is the "learning curve": a vector of the log likelihood values at each iteration.
%
% There are several optional arguments, which should be passed in the following order
%   LEARN_HMM(DATA, PRIOR, TRANSMAT, OBSMAT, MAX_ITER, THRESH, VERBOSE)
% These have the following meanings
%   max_iter = max. num EM steps to take (default 10)
%   thresh = threshold for stopping EM (default 1e-4)
%   verbose = 0 to suppress the display of the log lik at each iteration (Default 1).
%
% If the transition matrix is non-stationary (e.g., as in a POMDP),
% then TRANSMAT should be a cell array, where T{a}(i,j) = Pr(Q(t+1)=j|Q(t)=i,A(t)=a).
% The last arg should specify the sequence of actions in the same form as DATA:
%   LEARN_HMM(DATA, PRIOR, TRANSMAT, OBSMAT, MAX_ITER, THRESH, VERBOSE, As)
% The action at time 1 is ignored.
%
% If you want to clamp some of the parameters at fixed values, set the corresponding adjustable
% argument to 0 (default: everything is adjustable)
%   LEARN_HMM(..., VERBOSE, As, ADJ_PRIOR, ADJ_TRANS, ADJ_OBS)
%
% To avoid 0s when estimating OBSMAT, specify a non-zero equivalent sample size (e.g., 0.01) for
% the Dirichlet prior: LEARN_HMM(..., ADJ_OBS, DIRICHLET)
%
% When there is a single sequence, the smoothed posteriors using the penultimate set of
% parameters are returned in GAMMA:
%   [LL, PRIOR, TRANSMAT, OBSMAT, GAMMA] = LEARN_HMM(...)
% This can be useful for online learning and decision making.


dirichlet = 0;

previous_loglik = -inf;
loglik = 0;
converged = 0;
num_iter = 1;
LL = [];

numex = length(data);

while (num_iter <= max_iter) & ~converged
    % E step
    [loglik, Abar, pibar, Bbar, Cbar, Dbar, gamma] = ...
        compute_ess(prior, Amat, Bmat, Cmat, Dmat, data, dirichlet);
    
    if verbose, fprintf(1, 'iteration %d, loglik = %f\n', num_iter, loglik); end
    num_iter =  num_iter + 1;
    
    % M step
    prior = normalise(pibar);
    if ~isempty(Abar)
        Amat = mk_stochastic(Abar);
    end
    Bmat = mk_stochastic(Bbar);
    Cmat = mk_stochastic(Cbar);
    Dmat = mk_stochastic(Dbar);
    
    converged = em_converged(loglik, previous_loglik, thresh);
    previous_loglik = loglik;
    LL = [LL loglik];
end


%%%%%%%%%%%

function [loglik, Abar, pibar, Bbar, Cbar, Dbar, gamma] = ...
    compute_ess(prior, Amat, Bmat, Cmat, Dmat, data, dirichlet)
%
% Compute the Expected Sufficient Statistics for a discrete Hidden Markov Model.
%
% Outputs:
% exp_num_trans(i,j) = sum_l sum_{t=2}^T Pr(X(t-1) = i, X(t) = j| Obs(l))
% exp_num_visits1(i) = sum_l Pr(X(1)=i | Obs(l))
% exp_num_emit(i,o) = sum_l sum_{t=1}^T Pr(X(t) = i, O(t)=o| Obs(l))
% where Obs(l) = O_1 .. O_T for sequence l.

numex = length(data);
[S O] = size(Bmat);
Abar = zeros(S,S);
Cbar = zeros(S,S);
pibar = zeros(S,1);
Bbar = dirichlet*ones(S,O);
Dbar = dirichlet*ones(S,O);
loglik = 0;
estimated_trans = 0;

for ex=1:numex  % for each sequence
    obs = data{ex};
    obsV = obs(1,:);  % v -- the observations for this sequence (1xn)
    obsW = obs(2,:);  % w -- the observations for this sequence (1xn)
    T = length(obs);  % T -- the length of the sequence
    obslikB = mk_dhmm_obs_lik(obsV, Bmat);  % b_j -- the observation likelihoods (performance), (M_P x n)
    obslikC = mk_dhmm_obs_lik(obsW, Cmat);  % c_i -- the attention->knowledge likelihoods (attention), (M_A x n)
    obslikD = mk_dhmm_obs_lik(obsW, Dmat);  % d_j -- the observation likelihoods (attention), (M_A x n)
    [gamma, xi, current_ll] = forwards_backwards_ktattn(prior, Amat, obslikB, obslikC, obslikD);
    loglik = loglik +  current_ll;
    
    estimated_trans = 1;
    Abar = Abar + sum(xi,3);  % eqn (40b)
    for o=1:O
        ndx = find(obsW(1:end-1)==o);
        xi_w = sum(xi(:,:,ndx),1);
        Cbar(:,o) = Cbar(:,o) + sum(xi_w,3)';
    end
    
    pibar = pibar + gamma(:,1);  % eqn (40a)
    
    for o=1:O
        ndx = find(obsV==o);
        if ~isempty(ndx)  % eqn (40c)
            Bbar(:,o) = Bbar(:,o) + sum(gamma(:, ndx), 2);
        end
        ndx = find(obsW==o);
        if ~isempty(ndx)
            Dbar(:,o) = Dbar(:,o) + sum(gamma(:, ndx), 2);
        end
    end
end

if ~estimated_trans
    Abar = [];
end
