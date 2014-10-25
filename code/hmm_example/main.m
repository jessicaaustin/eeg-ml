
%% Generate observations and plot results

n = 100;
x = generateObservations(n);

colors = {'go', 'bo', 'ro'};
    
figure; hold on
    plot(1, sum(x==1), colors{1});
    plot(2, sum(x==2), colors{2});
    plot(3, sum(x==3), colors{3});
    legend('rainy', 'cloudy', 'sunny');
    xlabel('state');
    ylabel('count');
    xlim([0.5 3.5]);
    ylim([0 n]);
    grid on;
    title('counts for each state');
    
figure; hold on;
    for i=1:n
        plot(i, 1, colors{x(i)});
    end
    xlabel('time');
    ylabel('state');
    title(sprintf('states over time\nblue: cloudy, green: rainy, red: sunny'));
    
%% Estimate parameters of HMM
    
