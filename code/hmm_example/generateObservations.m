function observations = generateObservations(n)
% Generate n observations, based on a latent model

%% Internal Model

% Initial state distributions
p_init = [0.1;  % rainy
          0.5;  % cloudy
          0.4]; % sunny

% State transition probabilities
a11 = 0.3;  % rainy to rainy
a12 = 0.6;  % rainy to cloudy
a13 = 0.1;  % rainy to sunny
a22 = 0.7;  % cloudy to cloudy
a21 = 0.1;  % cloudy to rainy
a23 = 0.2;  % cloudy to sunny
a33 = 0.5;  % sunny to sunny
a31 = 0.1;  % sunny to rainy
a32 = 0.4;  % sunny to cloudy
A = [a11 a12 a13;
     a21 a22 a23;
     a31 a32 a33];

% each row should sum to one
assert(all((sum(A,2)-1)<eps));


%% Generate observations

observations = zeros(n,1);
observations(1) = discreteSample(p_init,1);
for i=1:n-1
    current_state = observations(i);
    state_transition_probabilities = A(current_state,:);
    next_state = discreteSample(state_transition_probabilities,1);
    observations(i+1) = next_state;
end
 
end