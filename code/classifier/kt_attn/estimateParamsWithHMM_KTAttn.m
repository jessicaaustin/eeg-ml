clear;
close all;

addpath('../../HMM_mat');

load('subjects.mat');
N = length(subjectids);

recompute = false;

%% Setup

% initial estimates for KT params
% a complete guess for the learning rate, and assume that the forgetting
% rate is near zero
p0 = [0.5; 
      0.5];
p_learn_0 = 0.7;
p_forget_0 = 0.01;
A0 =[1-p_forget_0  p_forget_0
     1-p_learn_0   p_learn_0];
p_guess_0 = 0.4;
p_slip_0 = 0.2;
B0 = [p_guess_0   1-p_guess_0
      p_slip_0    1-p_slip_0];
p_learn_a_0 = 0.7;
p_dontlearn_a_0 = 0.8;
C0 = [p_dontlearn_a_0   1-p_dontlearn_a_0
      1-p_learn_a_0     p_learn_a_0];
p_guess_a_0 = 0.6;
p_slip_a_0 = 0.7;
D0 = [p_slip_a_0      1-p_slip_a_0
      1-p_guess_a_0   p_guess_a_0];
       
% EM params
max_iter = 500;
thresh = 1e-8;
verbose = 1;


%% Estimate Params

for i=1:N
    sid=subjectids{i};
    fprintf('%s:\n', char(sid));

    filename = char(strcat('subjects/', sid, '_KTAttn.mat'));
    
    if ~recompute && exist(filename, 'file')
        fprintf('already exists, continuing\n');
        continue;
    end
    
    seqs_filename = char(strcat('subjects/', sid, '_sequences.mat'));
    load(seqs_filename);
    
    O = length(sequences.accept);
    x = cell(O,1);
    for o=1:O
        ASRobservations = sequences.accept{o};
        EEGobservations = binaryAttention(sequences.attention{o});
        x{o} = [ASRobservations';
                EEGobservations'];
    end
    
    [LL, p_est, A_est, B_est, C_est, D_est] = learn_dhmm_ktattn(x, p0, A0, B0, C0, D0, max_iter, thresh, verbose);

    l0 = p_est(2)
    p_learn =  A_est(2,2)
    p_forget = A_est(1,2)
    p_guess = B_est(1,1)
    p_slip = B_est(2,1)
    p_learn_a = C_est(2,2)
    p_dontlearn_a = C_est(1,1)
    p_guess_a = D_est(2,2)
    p_slip_a = D_est(1,1)
    
    subjectid = sid;
    p = p_est;
    A = A_est;
    B = B_est;
    C = C_est;
    D = D_est;
    save(filename, 'subjectid', ...
        'LL', 'A', 'B', 'C', 'D', 'p', ...
        'l0', 'p_learn', 'p_forget', 'p_guess', 'p_slip', ...
        'p_learn_a', 'p_dontlearn_a', 'p_guess_a', 'p_slip_a');

end


%% Display Results

N = length(subjectids);
all_l0 = zeros(N,1);
all_p_learn = zeros(N,1);
all_p_forget = zeros(N,1);
all_p_guess = zeros(N,1);
all_p_slip = zeros(N,1);
all_p_learn_a = zeros(N,1);
all_p_dontlearn_a = zeros(N,1);
all_p_guess_a = zeros(N,1);
all_p_slip_a = zeros(N,1);

for i=1:N
    sid=subjectids{i};
    filename = char(strcat('subjects/', sid, '_KTAttn.mat'));
    load(filename);
    
    all_l0(i) = l0;
    all_p_learn(i) = p_learn;
    all_p_forget(i) = p_forget;
    all_p_guess(i) = p_guess;
    all_p_slip(i) = p_slip;
    all_p_learn_a(i) = p_learn_a;
    all_p_dontlearn_a(i) = p_dontlearn_a;
    all_p_guess_a(i) = p_guess_a;
    all_p_slip_a(i) = p_slip_a;
end

figure; hold on;
    plot(repmat(p0(2), N, 1), 'k:');
    plot(all_l0,'ro')
    errorbar(round(N/2), mean(all_l0), std(all_l0));
    title('l0')
    ylim([0 1]);
    
figure; 
    subplot(1,2,1); hold on;
        plot(repmat(p_learn_0, N, 1), 'k:');
        plot(all_p_learn,'ro')
        errorbar(round(N/2), mean(all_p_learn), std(all_p_learn));
        title('p\_learn')
        ylim([0 1]);
    subplot(1,2,2); hold on;
        plot(repmat(p_forget_0, N, 1), 'k:');
        plot(all_p_forget,'ro')
        errorbar(round(N/2), mean(all_p_forget), std(all_p_forget));
        title('p\_forget')
        ylim([0 1]);
    suptitle('A')
    
figure;
    subplot(1,2,1); hold on;
        plot(repmat(p_guess_0, N, 1), 'k:');
        plot(all_p_guess,'ro')
        errorbar(round(N/2), mean(all_p_guess), std(all_p_guess));
        title('p\_guess')
        ylim([0 1]);
    subplot(1,2,2); hold on;
        plot(repmat(p_slip_0, N, 1), 'k:');
        plot(all_p_slip,'ro')
        errorbar(round(N/2), mean(all_p_slip), std(all_p_slip));
        title('p\_slip')
        ylim([0 1]);
    suptitle('B')

figure;
    subplot(1,2,1); hold on;
        plot(repmat(p_learn_a_0, N, 1), 'k:');
        plot(all_p_learn_a,'ro')
        errorbar(round(N/2), mean(all_p_learn_a), std(all_p_learn_a));
        title('p\_learn\_a')
        ylim([0 1]);
    subplot(1,2,2); hold on;
        plot(repmat(p_dontlearn_a_0, N, 1), 'k:');
        plot(all_p_dontlearn_a,'ro')
        errorbar(round(N/2), mean(all_p_dontlearn_a), std(all_p_dontlearn_a));
        title('p\_dontlearn\_a')
        ylim([0 1]);
    suptitle('C')
    
figure;
    subplot(1,2,1); hold on;
        plot(repmat(p_guess_a_0, N, 1), 'k:');
        plot(all_p_guess_a,'ro')
        errorbar(round(N/2), mean(all_p_guess_a), std(all_p_guess_a));
        title('p\_guess\_a')
        ylim([0 1]);
    subplot(1,2,2); hold on;
        plot(repmat(p_slip_a_0, N, 1), 'k:');
        plot(all_p_slip_a,'ro')
        errorbar(round(N/2), mean(all_p_slip_a), std(all_p_slip_a));
        title('p\_slip\_a')
        ylim([0 1]);
    suptitle('D')
    
    
%% Notes
    
% testing ideas:
% take a student
% use leave-one-out across sequences
% during each test sequence:
%   1. use l0 to estimate 1st state, and g,s to estimate 
%      observation based on that state
%   2. compare to actual observation
%   3. take actual observation. use l,f to estimate next state.
%      use g,s to estimate observation
%   4. compare to actual observation
%   5. repeat 3-4 for remaining elements in sequence
  