
clear;
close all;

addpath('../../HMM_mat');
addpath('../../Barber_ML');
addpath('../../Barber_ML_ext');

%% Generate observations and plot results

n = 10;
x = [];

xcolors = {'rx', 'gx'};
scolors = {'yo', 'ko'};
    
figure; hold on;

for w=1:100
    
    [x_w, S_w, A_true, b_true, p_true] = generateObservations(n);
    x(:,w) = x_w;
    
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
A_guess =[.5  .5
           0   1];
b_guess = [.5  .5
           .2  .8];
p_guess = [0.9; 
           0.1];

fprintf('Estimating parameters using EM (Baum Welch) MATLAB implementation:\n');
[A_est_EM_matlab, b_est_EM_matlab] = hmmtrain(x', A_guess, b_guess);

fprintf('Estimating parameters using EM (Baum Welch) HMM Toolbox implementation:\n');
[~, p_est_EM_HMMt, A_est_EM_HMMt, b_est_EM_HMMt] = learn_dhmm(x', p_guess, A_guess, b_guess, 'max_iter', 50);

fprintf('Estimating parameters using EM (Baum Welch) Barber implementation:\n');
opts.maxit = 30;
opts.plotprogress = 1;
x_cell = mat2cell(x,size(x,1),ones(1,size(x,2)));
figure;
[A_est_EM_Barber, p_est_EM_Barber, b_est_EM_Barber, ~] = HMMem(x_cell, size(A_guess,2), size(b_guess,2), opts);

%% Plot results

figure;
subplot(2,2,1); bar(p_true, 'r'); title('true initial p');
subplot(2,2,3); bar(p_est_EM_HMMt); title('learned initial p (HMM Toolbox)');
subplot(2,2,4); bar(p_est_EM_Barber); title('learned initial p (Barber)');

figure;
subplot(2,2,1); imagesc(A_true); title('true transition');
subplot(2,2,2); imagesc(A_est_EM_matlab); title('learned transition (MATLAB hmmtrain)');
subplot(2,2,3); imagesc(A_est_EM_HMMt); title('learned transition (HMM Toolbox)');
subplot(2,2,4); imagesc(A_est_EM_Barber); title('learned transition (Barber)');

figure;
subplot(2,2,1); imagesc(b_true); title('true emission');
subplot(2,2,2); imagesc(b_est_EM_matlab); title('learned emission (MATLAB hmmtrain)');
subplot(2,2,3); imagesc(b_est_EM_HMMt); title('learned emission (HMM Toolbox)');
subplot(2,2,4); imagesc(b_est_EM_Barber); title('learned emission (Barber)');


%% Print out results

Methods = {'Actual Value', 'MATLAB hmmtrain', 'HMM Toolbox', 'Barber'};

l0 = [p_true(2); 
      0;
      p_est_EM_HMMt(2);
      p_est_EM_Barber(2)];
p_learn =  [A_true(1,2); 
            A_est_EM_matlab(1,2); 
            A_est_EM_HMMt(1,2); 
            A_est_EM_Barber(1,2)];
p_forget = [A_true(2,1);
            A_est_EM_matlab(2,1);
            A_est_EM_HMMt(2,1);
            A_est_EM_Barber(2,1)];
p_guess = [b_true(1,2);
           b_est_EM_matlab(1,2);
           b_est_EM_HMMt(1,2);
           b_est_EM_Barber(1,2)];
p_slip = [b_true(2,1);
          b_est_EM_matlab(2,1);
          b_est_EM_HMMt(2,1);
          b_est_EM_Barber(2,1)];

T = table(l0, p_learn, p_forget, p_guess, p_slip, ...
    'RowNames', Methods)