
clear;
close all;

addpath('../../HMM_mat');
addpath('../../HMM_mat_ext');

%% Generate observations and plot results

n = 10;
x = {};

xcolors = {'rx', 'gx'};
scolors = {'yo', 'ko'};
    
figure; hold on;

for w=1:100
    
    [x_w, S_w, A_true, b_true, p_true] = generateObservations(n);
    x{w} = x_w';
    
    for i=1:n
        plot(i, 1+.5*w, scolors{S_w(i)}, 'LineWidth',3);
        plot(i, 1+.5*w, xcolors{x_w(i)}, 'LineWidth',2);
    end
end

xlabel('time');
ylabel('state (o), observation (x)');
title(sprintf('state and observations over time'));  
    
%% Estimate parameters of HMM
    
% Estimate parameters with just the sequence as input, using
% the Baum-Welch algorithm
A0 =[.5  .5
           0   1];
b0 = [.5  .5
           .2  .8];
p0 = [0.9; 
           0.1];

fprintf('Estimating parameters using EM (Baum Welch) HMM Toolbox ext implementation:\n');
[LL, p_est_EM_HMMt, A_est_EM_HMMt, b_est_EM_HMMt] = learn_dhmm_iohmm(x, p0, A0, b0, 50, 1E-8);

%% Plot results

figure;
subplot(1,2,1); bar(p_true, 'r'); title('true initial p');
subplot(1,2,2); bar(p_est_EM_HMMt); title('learned initial p (HMM Toolbox)');

figure;
subplot(1,2,1); imagesc(A_true); title('true transition');
subplot(1,2,2); imagesc(A_est_EM_HMMt); title('learned transition (HMM Toolbox)');

figure;
subplot(1,2,1); imagesc(b_true); title('true emission');
subplot(1,2,2); imagesc(b_est_EM_HMMt); title('learned emission (HMM Toolbox)');


%% Print out results

Methods = {'Actual Value', 'HMM Toolbox'};

l0 = [p_true(2); 
      p_est_EM_HMMt(2)];
p_learn =  [A_true(1,2); 
            A_est_EM_HMMt(1,2)];
p_forget = [A_true(2,1);
            A_est_EM_HMMt(2,1)];
p_guess = [b_true(1,2);
           b_est_EM_HMMt(1,2)];
p_slip = [b_true(2,1);
          b_est_EM_HMMt(2,1)];

T = table(l0, p_learn, p_forget, p_guess, p_slip, ...
    'RowNames', Methods)