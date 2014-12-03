function [LL, prior, Amat, Bmat, Cmat, Dmat, gamma] = learn_dhmm_ktattn(data, prior, Amat, Bmat, Cmat, Dmat, max_iter, thresh, ...
						  verbose, act, adj_prior, adj_trans, adj_obs, dirichlet)
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

%learn_dhmm(data, prior, transmat, obsmat, max_iter, thresh, verbose, act, adj_prior, adj_trans, adj_obs, dirichlet)
if nargin < 5+2, max_iter = 10; end
if nargin < 6+2, thresh = 1e-4; end
if nargin < 7+2, verbose = 1; end
if nargin < 8+2
  act = [];
  A = 0;
else
  A = length(Amat);
end
if nargin < 9+2, adj_prior = 1; end
if nargin < 10+2, adj_trans = 1; end
if nargin < 11+2, adj_obs = 1; end
if nargin < 12+2, dirichlet = 0; end

previous_loglik = -inf;
loglik = 0;
converged = 0;
num_iter = 1;
LL = [];

if ~iscell(data)
  data = num2cell(data, 2); % each row gets its own cell
end
if ~isempty(act) & ~iscell(act)
  act = num2cell(act, 2);
end
numex = length(data);


while (num_iter <= max_iter) & ~converged
  % E step
  [loglik, Abar, pibar, Bbar, Cbar, Dbar, gamma] = ...
      compute_ess(prior, Amat, Bmat, Cmat, Dmat, data, act, dirichlet);

  if verbose, fprintf(1, 'iteration %d, loglik = %f\n', num_iter, loglik); end
  num_iter =  num_iter + 1;

  % M step
  if adj_prior
    prior = normalise(pibar);
  end
  if adj_trans & ~isempty(Abar)
    if isempty(act)
      Amat = mk_stochastic(Abar);
    else
      for a=1:A
	Amat{a} = mk_stochastic(Abar{a});
      end
    end
  end
  if adj_obs
    Bmat = mk_stochastic(Bbar);
  end
  Cmat = mk_stochastic(Cbar);
  Dmat = mk_stochastic(Dbar);
  
  converged = em_converged(loglik, previous_loglik, thresh);
  previous_loglik = loglik;
  LL = [LL loglik];
end


%%%%%%%%%%%

function [loglik, Abar, pibar, Bbar, Cbar, Dbar, gamma] = ...
    compute_ess(prior, Amat, Bmat, Cmat, Dmat, data, act, dirichlet)
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
if isempty(act) % true
  Abar = zeros(S,S);
  Cbar = zeros(S,S);
  A = 0;
else
  A = length(Amat);
  Abar = cell(1,A);
  for a=1:A
    Abar{a} = zeros(S,S);
  end
  error('should not get here')
end
pibar = zeros(S,1);
Bbar = dirichlet*ones(S,O);
Dbar = dirichlet*ones(S,O);
loglik = 0;
estimated_trans = 0;

for ex=1:numex  % for each sequence
  obs = data{ex};
  obsV = obs(1,:);  % v -- the observations for this sequence (1xn)
  obsW = obs(2,:);  % v -- the observations for this sequence (1xn)
  T = length(obs);  % T -- the length of the sequence
  obslikB = mk_dhmm_obs_lik(obsV, Bmat);  % b_j -- the observation likelihoods (performance), (M_P x n)
  obslikC = mk_dhmm_obs_lik(obsW, Cmat);  % c_i -- the attention->knowledge likelihoods (attention), (M_A x n)
  obslikD = mk_dhmm_obs_lik(obsW, Dmat);  % d_j -- the observation likelihoods (attention), (M_A x n)
  if isempty(act)
    [gamma, xi, current_ll] = forwards_backwards_ktattn(prior, Amat, obslikB, obslikC, obslikD);
  else
    [gamma, xi, current_ll] = forwards_backwards_pomdp(prior, Amat, obslikB, act{ex});
  end
  loglik = loglik +  current_ll; 

  if T > 1
    estimated_trans = 1;
    if isempty(act) % true
      Abar = Abar + sum(xi,3);
    else
      % act(2) determines Q(2), xi(:,:,1) holds P(Q(1), Q(2))
      A = length(Amat);
      for a=1:A
        ndx = find(act{ex}(2:end)==a);
        if ~isempty(ndx)   % eqn (40b)
          Abar{a} = Abar{a} + sum(xi(:,:,ndx), 3);
        end
      end
      error('should not get here')
    end
    
    for o=1:O
        ndx = find(obsW(1:end-1)==o);
        xi_w = sum(xi(:,:,ndx),2);
        Cbar(:,o) = Cbar(:,o) + sum(xi_w,3);
    end
  end
              % eqn (40a)
  pibar = pibar + gamma(:,1);
  
  if T < O  % false
    for t=1:T
      o = obsV(t);
      Bbar(:,o) = Bbar(:,o) + gamma(:,t);
    end
    error('should not get here (sequence less than num states)');
  else
    for o=1:O
      ndx = find(obsV==o);
      if ~isempty(ndx)  % eqn (40c)
        Bbar(:,o) = Bbar(:,o) + sum(gamma(:, ndx), 2);
      end
      ndx = find(obsW==o);
      if ~isempty(ndx)  % eqn (40c)
        Dbar(:,o) = Dbar(:,o) + sum(gamma(:, ndx), 2);
      end
    end
  end
end

if ~estimated_trans
  Abar = [];
end
