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
% have to rely on observed data. Currently, we rely on the barometric
% pressure to get an idea of the current weather.
% The barometer has 3 possible readings: low, medium, and high pressure.
% Generally, lower pressure means rain, medium pressure means
% cloudy, and high pressure means sun.
% The probabilities of low, medium, and high readings for each
% weather state is given by the matrix b.
%


% Initial state distributions
p = [0.8;  % unknown
     0.2];  % known

% State transition probabilities
a11 = 0.3;  % unknown to unknown (they didn't learn anything this step)
a12 = 0.7;  % unknown to known   (they learned something!)
a22 = 1.0;  % known to known   (we assume forget_rate=0)
a21 = 0.0;  % known to unknown (we assume forget_rate=0)
A = [a11 a12;
     a21 a22]

% each row should sum to one
assert(all((sum(A,2)-1)<eps));

% observation probability distribution for each state
   % reject  accept   <--ASR result
b = [.6       .4;  % unknown
     .1       .9]  % known

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