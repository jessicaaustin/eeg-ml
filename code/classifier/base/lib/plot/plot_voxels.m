function h = plot_voxels(meta, info, data)

col = mean(cell2mat(data));
coord = NaN * ones(size(meta.coordToCol));
for c = 1:length(col),
	p = meta.colToCoord(c,:);
	coord(p(1), p(2), p(3)) = col(c);
end

h = gcf;
colorbar;

c = 0;
for plane = 1:meta.dimz,
	c = c + 1;
	subplot(5,5, c);

	X = coord(:,:,plane);

	h2 = imagesc(X, [-1 1]);
	set(gca, 'XTickLabel', {});
	set(gca, 'YTickLabel', {});
end
