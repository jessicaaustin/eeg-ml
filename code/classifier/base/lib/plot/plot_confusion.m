function h = plot_confusion(scores, Y, labels);

% renormalizing the scores to make the plot more revealing
scores = scores - repmat(min(scores,[],2), 1, size(scores,2)); 
scores = scores ./ repmat(sum(scores,2), 1, size(scores,2));

[tY I] = sort(Y);
scores = scores(I,:);
labels = labels(I);

h = imagesc(scores);
colorbar;
set(gca, 'XTick',      1:length(labels));
set(gca, 'XTickLabel', labels);
set(gca, 'YTick',      1:length(labels));
set(gca, 'YTickLabel', labels);

angle = 90;
rotateticklabel(gca,angle);
