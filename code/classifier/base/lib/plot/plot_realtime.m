% the data
     nt = 5; % <- traces
     np = 100; % <- data/trace

% prepare the plot
     axes('xlim', [1, np], 'ylim', [-2, 5]);
     x = 1:np;
     y = -inf * ones(size(x));
     lh = line(x, y ,...
            'marker', '.', ...
            'markersize', 5, ...
            'linestyle', 'none');
     lb = line([inf, inf], [-2, 5]);
     shg;

% gather the data and plot in <real-time>...
for i = 1:nt*np
     ix = rem(i - 1, np) + 1;
     y(ix) = .5 * fix(i / np) + rand; % <- new data
     set(lh, 'ydata', y);
     set(lb, 'xdata', [ix, ix]);
     pause(.0001); % <- a time consuming OP
end
