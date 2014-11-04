
clear;
close all;

code_root = ('../../');
addpath([code_root, 'common']);

%% Generate observations and plot results

n = 1000;
[x, S] = generateObservations(n);

colors = {'co', 'bo', 'ko'};
    
figure; hold on
    plot(1, mean(x(S==1)), colors{1});
    plot(2, mean(x(S==2)), colors{2});
    plot(3, mean(x(S==3)), colors{3});
    xlabel('state');
    ylabel('mean pressure reading');
    xlim([0.5 3.5]);
    ylim(ylim.*[.5 1.1])
    grid on;
    title('mean pressure reading for each state');
    
%% Estimate parameters of HMM
    
% Estimate parameters using training data (with a known sequence),
% using MLE. Note: we're assuming that the observations are conditionally
% independent given a state, and thus estimate them separately.
% (See this SO post for an idea of what I mean:
% http://stackoverflow.com/questions/17487356/hidden-markov-model-for-multiple-observed-variables)
fprintf('Estimating parameters using MLE and training data with known states, MATLAB implementation:\n');
A_est_MLE_all = [];
b_est_MLE = zeros(3,3,size(x,2));
for i=1:size(x,2)  % iterate over all the observations
    [A_est_i, b_est_i] = hmmestimate(x(:,i), S);
    A_est_MLE_all(:,:,i) = A_est_i;
    b_est_MLE(:,:,i) = b_est_i;
end
A_est_MLE = mean(A_est_MLE_all,3)
b_est_MLE

