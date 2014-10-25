
clear;
close all;

%% Generate observations and plot results

n = 10000;
[x, S, mx, mS] = generateObservations(n);

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
    
