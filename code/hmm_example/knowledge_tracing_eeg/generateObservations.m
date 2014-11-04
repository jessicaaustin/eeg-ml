function [observations, hiddenStates, A, b, p] = generateObservations(n)
% Generate n observations, based on a latent model

%% Internal Model
%
% A particular word has two states:
%  S1: known
%  S2: unknown
% The initial state distribution is given as p_init, and 
% the state transition probabilities are given in A.
%
% Of course, we can't observe these states directly, but instead
% have to rely on observed data. Currently, we rely on the ASR to give
% us an idea of whether or not the student got the word correct. For
% each task, the ASR will return either "incorrect" or "correct". 
%


% Initial state distributions
l0 = 0.2;
p = [1-l0;  % unknown
     l0];  % known

% State transition probabilities
p_learn = 0.7;
p_forget = 0.05;
A = [1-p_learn   p_learn;
     p_forget    1-p_forget];

% each row should sum to one
assert(all((sum(A,2)-1)<eps));

% observation probability distribution for each state
p_guess = .4;
p_slip = .2;
   % incorrect  correct  <-- observed result
b = [1-p_guess  p_guess;   % unknown
     p_slip    1-p_slip];  % known
                          %   ^ knowledge state 

%% Generate observations

% generate states
states = zeros(n,1);
states(1) = discreteSample(p,1);
for i=1:n-1
    current_state = states(i);
    state_transition_probabilities = A(current_state,:);
    next_state = discreteSample(state_transition_probabilities,1);
    states(i+1) = next_state;
end

% generate observations for each state
observations = zeros(n,1);
for i=1:n
    sensor_probabilities_for_state = b(states(i),:);
    observations(i) = discreteSample(sensor_probabilities_for_state,1);
end
 
%% Return state
% Of course, this is a latent variable so we couldn't directly 
% know it. However, return this variable for debugging purposes.

hiddenStates = states;

end