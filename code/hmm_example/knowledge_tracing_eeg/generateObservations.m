function [ASRobservations, EEGobservations, attentionStates, knowledgeStates, ...
          p0_KT, A_KT, b_asr] = generateObservations(n)
% Generate n observations, based on a latent model

%% Internal Model: Knowledge State

% A particular word has two states:
%  S1: unknown
%  S2: known
% The initial state distribution is given as p_init, and 
% the state transition probabilities are given in A.
%

% Initial state distributions
l0 = 0.2;
p0_KT = [1-l0;  % unknown
         l0];  % known

% State transition probabilities
p_learn = 0.7;
p_forget = 0.05;
      %   un->un     un->k
A_KT = [1-p_learn   p_learn;
      %   k->un     k->k
        p_forget    1-p_forget];

% each row should sum to one
assert(all((sum(A_KT,2)-1)<eps));

% Of course, we can't observe these states directly, but instead
% have to rely on observed data. Currently, we rely on the ASR to give
% us an idea of whether or not the student got the word correct. For
% each task, the ASR will return either "incorrect" or "correct". 
%

% observation probability distribution for each state
p_guess = .4;
p_slip = .2;
   % incorrect  correct  <-- observed result
b_asr = [1-p_guess  p_guess;   % unknown
         p_slip    1-p_slip];  % known
                             %   ^ knowledge state 

%% Internal Model: Attention State
%
% At any particular time, a student could be in one of two states:
%  S1: inattentive
%  S2: attentive
% The probability of being in either state is influenced by the time
% elapsed since the last task, T. Specifically, there is some optimal time
% Tbest that we wish to learn. Any deviation from Tbest decreases the
% probability that the student is attentive.
%

% p(attentive) = N(T_best, s_attention)
T_best = 1; %[second]
s_attention = 0.1;  %[variation]

% Knowledge state transitions based on attention state
p_ak = .8;  % probability that if you pay attention then you learn something
         % unknown  known
A_atten = [p_ak   1-p_ak;   % inattentive
           1-p_ak   p_ak];      % attentive

% We observe the time that was elapsed since the last task.
% We also observe EEG signal, which gives us an observation of the
% attention state.
eeg_false_pos = 0.2;
eeg_false_neg = 0.3;
            % inattentive  attentive  <-- eeg signal
b_eeg = [1-eeg_false_pos  eeg_false_pos;    % inattentive
         eeg_false_neg   1-eeg_false_neg];  % attentive
                                           %   ^ attention state 
 
ac = .2;

%% Generate Attention states and observations

% randomly generate times between questions
breakTimes = abs(randn(n,1) + 1.5);

% probability that the student is attentive, based on the break time
p_attentive = normalPDF(T_best, s_attention, breakTimes);
% threshold to get states
attentionStates = p_attentive > .5;
% state is either 1 or 2
attentionStates = attentionStates + 1;

% generate EEG signals
EEGobservations = zeros(n,1);
for i=1:n
    eeg_probabilities_for_state = b_eeg(attentionStates(i),:);
    EEGobservations(i) = discreteSample(eeg_probabilities_for_state,1);
end

%% Generate Knowledge states and observations

% generate states
knowledgeStates = zeros(n,1);
knowledgeStates(1) = discreteSample(p0_KT,1);
for i=1:n-1
    current_state = knowledgeStates(i);
    current_attention_state = attentionStates(i);
    state_transition_probabilities = A_KT(current_state,:);
    attention_probabilities = A_atten(current_attention_state,:);
    total_transition_prob = sum(state_transition_probabilities'*attention_probabilities,2);
    next_state = discreteSample(total_transition_prob,1);
    knowledgeStates(i+1) = next_state;
end

% generate observations for each state
ASRobservations = zeros(n,1);
for i=1:n
    sensor_probabilities_for_state = b_asr(knowledgeStates(i),:);
    ASRobservations(i) = discreteSample(sensor_probabilities_for_state,1);
end
 
end