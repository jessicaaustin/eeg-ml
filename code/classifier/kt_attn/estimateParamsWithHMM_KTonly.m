clear;
close all;

addpath('../../HMM_mat');

load('subjects.mat');
N = length(subjectids);

recompute = false;

%% Estimate Params

for i=1:N
    subjectid=subjectids{i};
    fprintf('%s:\n', char(subjectid));

    filename = char(strcat('subjects/', subjectid, '_KTparams.mat'));
    
    if ~recompute && exist(filename, 'file')
        fprintf('already exists, continuing\n');
        continue;
    end
    
    seqs_filename = char(strcat('subjects/', subjectid, '_sequences.mat'));
    load(seqs_filename);
    
    [LL,p,A,B,C,D, l0, p_learn, p_forget, p_guess, p_slip] ...
        = estimateParamsForSubject(sequences, 1:length(sequences), 'KT');
   
    save(filename, 'subjectid', ...
        'LL', 'A', 'B', ...
        'l0', 'p_learn', 'p_forget', 'p_guess', 'p_slip');

end


%% Display Results

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
  
N = length(subjectids);
all_l0 = zeros(N,1);
all_p_learn_KT = zeros(N,1);
all_p_forget_KT = zeros(N,1);
all_p_guess_KT = zeros(N,1);
all_p_slip_KT = zeros(N,1);

for i=1:N
    subjectid=subjectids{i};
    filename = char(strcat('subjects/', subjectid, '_KTparams.mat'));
    load(filename);
    
    all_l0(i) = l0;
    all_p_learn_KT(i) = p_learn;
    all_p_forget_KT(i) = p_forget;
    all_p_guess_KT(i) = p_guess;
    all_p_slip_KT(i) = p_slip;
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
        plot(all_p_learn_KT,'ro')
%         errorbar(round(N/2), mean(all_p_learn), std(all_p_learn));
        title('Learn Rate (l)')
        ylim([0 1]);
    subplot(1,2,2); hold on;
        plot(repmat(p_forget_0, N, 1), 'k:');
        plot(all_p_forget_KT,'ro')
%         errorbar(round(N/2), mean(all_p_forget), std(all_p_forget));
        title('Forget Rate (f)')
        ylim([0 1]);
    suptitle('A')
    
figure;
    subplot(1,2,1); hold on;
        plot(repmat(p_guess_0, N, 1), 'k:');
        plot(all_p_guess_KT,'ro')
%         errorbar(round(N/2), mean(all_p_guess), std(all_p_guess));
        title('p\_guess')
        ylim([0 1]);
    subplot(1,2,2); hold on;
        plot(repmat(p_slip_0, N, 1), 'k:');
        plot(all_p_slip_KT,'ro')
%         errorbar(round(N/2), mean(all_p_slip), std(all_p_slip));
        title('p\_slip')
        ylim([0 1]);
    suptitle('B')


    
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
  