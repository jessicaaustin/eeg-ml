
clear;
close all;

addpath('../../common');
addpath('../../HMM_mat');
addpath('../../HMM_mat_ext');

verbose = 1;
maxiter = 50;
thresh = 1e-8;

%% Generate observations

n = 50;  % number of times they saw the same word
x = {};

for w=1:100  % number of words they saw
    
    [ASRobservations, EEGobservations, attentionStates, knowledgeStates, ...
        l0_true, p_learn_true, p_forget_true, p_guess_true, p_slip_true,  ...
        p_learn_a_true, p_dontlearn_a_true, p_guess_a_true, p_slip_a_true] = generateObservations(n);
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
d0 = [.6  .4
      .4  .6];

fprintf('Estimating parameters using EM (Baum Welch) HMM Toolbox modified implementation:\n');
[LL, p_est, A_est, B_est, C_est, D_est] = learn_dhmm_ktattn(x, p0, A0, b0, c0, d0, maxiter, thresh, verbose);


%% Print out results

Methods = {'Actual Value', 'HMM Toolbox'};

l0 = [l0_true; 
      p_est(2)];
p_learn =  [p_learn_true; 
            A_est(2,2)];
p_forget = [p_forget_true;
            A_est(1,2)];
p_guess = [p_guess_true;
           B_est(1,1)];
p_slip = [p_slip_true;
          B_est(2,1)];
p_learn_a = [p_learn_a_true;
             C_est(2,2)];
p_dontlearn_a = [p_dontlearn_a_true;
                 C_est(1,1)];
p_guess_a = [p_guess_a_true;
             D_est(2,2)];
p_slip_a = [p_slip_a_true;
             D_est(1,1)];

T = table(l0, p_learn, p_forget, p_guess, p_slip, p_learn_a, p_dontlearn_a, p_guess_a, p_slip_a, ...
    'RowNames', Methods)