
clear;
close all;

code_root = ('../../');
addpath([code_root, 'common']);
addpath([code_root, 'HMM_mat']);
addpath([code_root, 'Barber_ML']);
addpath([code_root, 'Barber_ML_ext']);

%% Generate observations and plot results

n = 1000;
[x, S, mx, mS, A_true, b_true, p_true] = generateObservations(n);

colors = {'co', 'bo', 'ko'};
mcolors = {'cx', 'bx', 'kx'};
    
figure; hold on
    plot(1, mean(x(S==1)), colors{1});
    plot(2, mean(x(S==2)), colors{2});
    plot(3, mean(x(S==3)), colors{3});
    plot(1, mean(mx(mS==1)), mcolors{1});
    plot(2, mean(mx(mS==2)), mcolors{2});
    plot(3, mean(mx(mS==3)), mcolors{3});
    xlabel('state');
    ylabel('mean pressure reading');
    xlim([0.5 3.5]);
    ylim(ylim.*[.5 1.1])
    grid on;
    title('mean pressure reading for each state');
    
if n<=100    
    figure; hold on;
        for i=1:n
            plot(i, 1, colors{x(i)});
            plot(i, 1.1, mcolors{mx(i)});
            ylim([.5 1.5]);
        end
        xlabel('time');
        ylabel('reading');
        title(sprintf('readings over time'));  
end
    
%% Estimate parameters of HMM
    
% Estimate parameters with just the sequence as input, using
% the Baum-Welch algorithm
A_guess = [.50 .25 .25;
           .25 .50 .25;
           .25 .25 .50];
b_guess = [.8 .1 .1;
           .1 .8 .1;
           .1 .1 .8];
p_guess = [1/3 1/3 1/3];


fprintf('Estimating parameters using EM (Baum Welch) MATLAB implementation:\n');
[A_est_EM_matlab, b_est_EM_matlab] = hmmtrain(x, A_guess, b_guess)

fprintf('Estimating parameters using EM (Baum Welch) HMM Toolbox implementation:\n');
[~, p_est_EM_HMMt, A_est_EM_HMMt, b_est_EM_HMMt] = learn_dhmm(x, p_guess, A_guess, b_guess, 'max_iter', 50)

fprintf('Estimating parameters using EM (Baum Welch) Barber implementation:\n');
opts.maxit = 30;
opts.plotprogress = 1;
x_cell{1} = x;
figure;
[A_est_EM_Barber, p_est_EM_Barber, b_est_EM_Barber, ~] = HMMem(x_cell, 3, 3, opts)

fprintf('Estimating parameters using EM (Baum Welch) Barber implementation with initial guess:\n');
figure;
[A_est_EM_BarberWithGuess, p_est_EM_BarberWithGuess, b_est_EM_BarberWithGuess, ~] = HMMem_withInitGuess(x_cell, A_guess, b_guess, p_guess', opts)

% Estimate parameters using training data (with a known sequence),
% using MLE

fprintf('Estimating parameters using MLE and training data with known states, MATLAB implementation:\n');
[A_est_MLE, b_est_MLE] = hmmestimate(x, S)

%% Plot results

figure;
subplot(2,2,1); bar(p_true, 'r'); title('true initial p');
subplot(2,2,2); bar(p_est_EM_HMMt); title('learned initial p (HMM Toolbox)');
subplot(2,2,3); bar(p_est_EM_HMMt); title('learned initial p (Barber)');
subplot(2,2,4); bar(p_est_EM_HMMt); title('learned initial p (Barber with guess)');

figure;
subplot(2,3,1); imagesc(A_true); title('true transition');
subplot(2,3,2); imagesc(A_est_EM_matlab); title('learned transition (MATLAB hmmtrain)');
subplot(2,3,3); imagesc(A_est_EM_HMMt); title('learned transition (HMM Toolbox)');
subplot(2,3,4); imagesc(A_est_EM_Barber); title('learned transition (Barber)');
subplot(2,3,5); imagesc(A_est_EM_BarberWithGuess); title('learned transition (Barber with guess)');
subplot(2,3,6); imagesc(A_est_MLE); title('learned transition (MLE)');

figure;
subplot(2,3,1); imagesc(b_true); title('true emission');
subplot(2,3,2); imagesc(b_est_EM_matlab); title('learned emission (MATLAB hmmtrain)');
subplot(2,3,3); imagesc(b_est_EM_HMMt); title('learned emission (HMM Toolbox)');
subplot(2,3,4); imagesc(b_est_EM_Barber); title('learned emission (Barber)');
subplot(2,3,5); imagesc(b_est_EM_BarberWithGuess); title('learned emission (Barber with guess)');
subplot(2,3,6); imagesc(b_est_MLE); title('learned emission (MLE)');


