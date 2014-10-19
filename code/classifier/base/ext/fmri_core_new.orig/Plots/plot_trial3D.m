%
% In:
% - a matrix containing one image to plot per row (could be the
%   matrix for a trial in the data cell array)
% - the meta for that data
%
% Out:
% - interactive plot in 3D of each image
%
% Notes:
% - this function is a skeleton that should be easy to modify,
%   hence it's not printing any information about subject, study, etc
%
% Examples:
% - plot_trial3D(data{1},meta)
% - plot_trial3D(imageToPlot,meta)
%
% History
% - 2005 Jul 23 - fpereira - created
%

function [] = plot_trial3D(varargin)

%
% Process parameters and setup
%

this = 'plot_trial3D';
l = length(varargin);
if l < 2; eval(sprintf('help %s',this)); return; else

  matrixToPlot = varargin{1};
  meta         = varargin{2};
  
  if l > 2; scale = varargin{3}; else
    % compute an appropriate scale for the data
    scale = prctile(matrixToPlot(:),[1 99]);
  end
end

% find indices in 3D for all voxels
dimensions  = [meta.dimx,meta.dimy,meta.dimz];
c = meta.colToCoord;
indicesIn3D = sub2ind(dimensions,c(:,1),c(:,2),c(:,3));

% each image will be placed in this volume
volume = repmat(NaN,dimensions);
nrows = ceil(sqrt(meta.dimz)); ncols = nrows;

%
% Plot the images
%

for t = 1:size(matrixToPlot,1)

  fprintf('plotting image %d\n',t);
  volume(indicesIn3D) = matrixToPlot(t,:);

  for z = 1:meta.dimz
    
    subplot(nrows,ncols,z);
    
    imagesc(rot90(volume(:,:,z),-1),scale);
    if z == 1; colorbar('vertical'); end
    view(0,90);
  
    xlabel(sprintf('z=%d',z));
    set(gca,'XTick',[]);set(gca,'YTick',[]);axis square;
  end
  pause
end
