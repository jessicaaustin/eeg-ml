
clear;
close all;

addpath('../../common');
addpath('../../HMM_mat');
addpath('../../HMM_mat_ext');

%% Generate observations and plot results

n = 10;
[ASRobservations, EEGobservations, attentionStates, knowledgeStates, l0_true, p_learn_true, ...
    p_forget_true, p_guess_true, p_slip_true, p_ak_true, p_ap_true] = generateObservations(n);

%% Estimate parameters of HMM
    
% Estimate parameters with just the sequence as input, using
% the Baum-Welch algorithm
A0 =[.5  .5
           0   1];
b0 = [.5  .5
           .2  .8];
p0 = [0.9; 
           0.1];

x = ASRobservations;
fprintf('Estimating parameters using EM (Baum Welch) HMM Toolbox ext implementation:\n');
[LL, p_est_EM_HMMt, A_est_EM_HMMt, b_est_EM_HMMt] = learn_dhmm_iohmm(x, p0, A0, b0, 50, 1E-8);


%% Print out results

Methods = {'Actual Value', 'HMM Toolbox'};

l0 = [l0_true; 
      p_est_EM_HMMt(2)];
p_learn =  [p_learn_true; 
            A_est_EM_HMMt(1,2)];
p_forget = [p_forget_true;
            A_est_EM_HMMt(2,1)];
p_guess = [p_guess_true;
           b_est_EM_HMMt(1,2)];
p_slip = [p_slip_true;
          b_est_EM_HMMt(2,1)];
p_ak = [p_ak_true;
        NaN];
p_ap = [p_ap_true;
        NaN];

T = table(l0, p_learn, p_forget, p_guess, p_slip, p_ak, p_ap, ...
    'RowNames', Methods)