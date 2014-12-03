
clear;
close all;

addpath('../../common');
addpath('../../HMM_mat');
addpath('../../HMM_mat_ext');

%% Generate observations

n = 10;  % number of times they saw the same word
x = {};

for w=1:100  % number of words they saw
    
    [ASRobservations, EEGobservations, attentionStates, knowledgeStates, ...
        l0_true, p_learn_true, p_forget_true, p_guess_true, p_slip_true, p_learn_a_true, p_guess_a_true] = generateObservations(n);
    x{w} = [ASRobservations';
            EEGobservations'];
    
end

%% Estimate parameters of HMM
    
% Estimate parameters with just the sequence as input, using
% the Baum-Welch algorithm
p0 = [0.9; 
      0.1];
A0 =[.5  .5
    0.01   0.99];
b0 = [.5  .5
      .2  .8];
c0 = [.6  .4
      .4  .6];
d0 = [.6  .6
      .4  .4];

fprintf('Estimating parameters using EM (Baum Welch) HMM Toolbox ext implementation:\n');
[LL, p_est, A_est, B_est, C_est, D_est] = learn_dhmm_ktattn(x, p0, A0, b0, c0, d0, 50, 1E-8, 0);


%% Print out results

Methods = {'Actual Value', 'HMM Toolbox'};

l0 = [l0_true; 
      p_est(2)];
p_learn =  [p_learn_true; 
            A_est(1,2)];
p_forget = [p_forget_true;
            A_est(2,1)];
p_guess = [p_guess_true;
           B_est(1,2)];
p_slip = [p_slip_true;
          B_est(2,1)];
p_ak = [p_learn_a_true;
        C_est(2,1)];
p_ap = [p_guess_a_true;
        D_est(2,1)];

T = table(l0, p_learn, p_forget, p_guess, p_slip, p_ak, p_ap, ...
    'RowNames', Methods)