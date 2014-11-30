clear;

addpath('../../HMM_mat');

load('hmmData_2012_2013');
subject = hmmData_2012_2013{1};

x = subject.fluent;
for i=1:length(x)
    x{i} = x{i} + ones(size(x{i}));
end

A0 =[.5  .5
           0   1];
b0 = [.5  .5
           .2  .8];
p0 = [0.5; 
           0.5];

fprintf('Estimating parameters using EM (Baum Welch) HMM Toolbox implementation:\n');
[LL, p_est_EM_HMMt, A_est_EM_HMMt, b_est_EM_HMMt] = learn_dhmm(x, p0, A0, b0, 50, 1E-8);

l0 = p_est_EM_HMMt(2)
p_learn =  A_est_EM_HMMt(1,2)
p_forget = A_est_EM_HMMt(2,1)
p_guess = b_est_EM_HMMt(1,2)
p_slip = b_est_EM_HMMt(2,1)


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
  