clear;
close all;

addpath('../../HMM_mat');

load('subjects.mat');

%% Estimate Params

% initial estimates for KT params
% a complete guess for the learning rate, and assume that the forgetting
% rate is near zero
A0 =[.5  .5
     0.01   0.99];
b0 = [.5  .5
      .5  .5];
p0 = [0.5; 
           0.5];
       
% EM params
max_iter = 500;
thresh = 1e-8;

for sid=subjectids'
    fprintf('%s:\n', char(sid));

    filename = char(strcat('subjects/', sid, '_KTparams.mat'));
    
    if exist(filename, 'file')
        fprintf('already exists, continuing\n');
        continue;
    end
    
    seqs_filename = char(strcat('subjects/', sid, '_sequences.mat'));
    load(seqs_filename);
    
    x = sequences.accept;
    
    [LL, p_est, A_est, b_est] = learn_dhmm(x, p0, A0, b0, max_iter, thresh);

    l0 = p_est(2)
    p_learn =  A_est(1,2)
    p_forget = A_est(2,1)
    p_guess = b_est(1,2)
    p_slip = b_est(2,1)

    subjectid = sid;
    A = A_est;
    b = b_est;
    save(filename, 'subjectid', ...
        'LL', 'A', 'b', ...
        'l0', 'p_learn', 'p_forget', 'p_guess', 'p_slip');

end


%% Display Results

N = length(subjectids);
all_l0 = zeros(N,1);
all_p_learn = zeros(N,1);
all_p_forget = zeros(N,1);
all_p_guess = zeros(N,1);
all_p_slip = zeros(N,1);

for i=1:N
    sid=subjectids{i};
    filename = char(strcat('subjects/', sid, '_KTparams.mat'));
    load(filename);
    
    all_l0(i) = l0;
    all_p_learn(i) = p_learn;
    all_p_forget(i) = p_forget;
    all_p_guess(i) = p_guess;
    all_p_slip(i) = p_slip;
end

figure; hold on;
    plot(repmat(p0(2), N, 1), 'k:');
    plot(all_l0,'ro')
    errorbar(round(N/2), mean(all_l0), std(all_l0));
    title('l0')
    ylim([0 1]);
figure; hold on;
    plot(repmat(A0(1,2), N, 1), 'k:');
    plot(all_p_learn,'ro')
    errorbar(round(N/2), mean(all_p_learn), std(all_p_learn));
    title('p\_learn')
    ylim([0 1]);
figure; hold on;
    plot(repmat(A0(2,1), N, 1), 'k:');
    plot(all_p_forget,'ro')
    errorbar(round(N/2), mean(all_p_forget), std(all_p_forget));
    title('p\_forget')
    ylim([0 1]);
figure; hold on;
    plot(repmat(b0(1,2), N, 1), 'k:');
    plot(all_p_guess,'ro')
    errorbar(round(N/2), mean(all_p_guess), std(all_p_guess));
    title('p\_guess')
    ylim([0 1]);
figure; hold on;
    plot(repmat(b0(2,1), N, 1), 'k:');
    plot(all_p_slip,'ro')
    errorbar(round(N/2), mean(all_p_slip), std(all_p_slip));
    title('p\_slip')
    ylim([0 1]);

    
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
  