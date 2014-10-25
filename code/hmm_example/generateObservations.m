function [observations, hiddenStates, mx, mS] = generateObservations(n)
% Generate n observations, based on a latent model

%% Internal Model
%
% This weather model has three states
%  S1: rainy
%  S2: cloudy
%  S3: sunny
% The initial state distribution is given as p_init, and 
% the state transition probabilities are given in A.
%
% Of course, we can't observe these states directly, but instead
% have to rely on sensor data. Currently, we rely on the barometric
% pressure to get an idea of the current weather.
% The barometer has 3 possible readings: low, medium, and high pressure.
% Generally, lower pressure means rain, medium pressure means
% cloudy, and high pressure means sun.
% The probabilities of low, medium, and high readings for each
% weather state is given by the matrix b.
%

% TODO: add more sensors:
%  v1: air pressure (barometer)
%  v2: temperature (thermometer)
%  v3: relative humidity (hygrometer)
%  v4: precipitation (rain gauge)
%  v5: wind speed (anemometer)
% TODO: make sensors noisy: draw from a gaussian


% Initial state distributions
p = [0.1;  % rainy
     0.5;  % cloudy
     0.4]  % sunny

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
     a31 a32 a33]

% each row should sum to one
assert(all((sum(A,2)-1)<eps));

% observation probability distribution for each state
   % low  medium  high pressure reading
b = [.7    .2      .1;  % rainy
     .3    .5      .2;  % cloudy
     .1    .4      .5] % sunny

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

%% Sample using the HMM toolbox implementation
% for comparison purposes. note that there is no way to specific the
% initial state (you can specify initial state distribution though), so
% this will necessarily deviate from the above result (x and S)
[mx, mS] = sample_dhmm(p, A, b, 1, n);  

%% Sample using the MATLAB implementation
% [mx,mS] = hmmgenerate(n,A,b);


end