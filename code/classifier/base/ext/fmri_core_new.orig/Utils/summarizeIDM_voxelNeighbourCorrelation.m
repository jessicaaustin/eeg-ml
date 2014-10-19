% Computes maps of average correlation of a voxel with its neighbours, during each condition.
% In more detail, for each condition:
% - take all trials in this condition and concatenate them,
%   yielding one time series per voxel.
% - for each voxel, compute the average correlation of that time
%   series with the time series of each voxel that borders it in 3D
%
% In:
% - info+data+meta
%
% Out:
% - Three arrays with #conditions x #voxels, where each row
%   pertains to one condition, in the order of sorted condition labels
%   - correlationImages -  each row is the image containing the
%     correlation values of each voxel with its neighbours.
%   - pvalueImages - each row is an image containing p-values for
%     the correlation values in the corresponding row of
%     correlationImages (this is a bit dodgy, as the p-value
%     for the correlation between two voxels, under the null
%     hypothesis that there is no correlation, is what is computed
%     for all the neighbours of a voxel. It's then averaged over
%     neighbours. To do this right we'd need the distribution of
%     the mean correlation with neighbours, assuming there's no
%     correlation with any of them. But I don't have that distribution...)
%   - meanImages - each row contains the image of mean values of
%     the time series of each voxel used to compute correlations
%   - stdvImages - same thing for standard deviations
%
% Examples:
%
% [correlationImages,pvalueImages,meanImages,stdvImages,conditions] = summarizeIDM_voxelNeighbourCorrelation(info,data,meta);
%  
% plot_trial3D( correlationImages(2,:), meta );
% (plots the correlation image for condition # <conditions(2)>)
%
%
% Notes:
% - Condition 0 is ignored
% - This can be made to use only specific trials by modifying the
%   info in the IDM passed in (e.g. zero out the condition #s of trials to exclude them)
% - If you want to plot the images one of the examples uses "plot_trial3D" to do it 
% - This function should be easily modifiable to do any computation over the neighbours of a voxel

%
% History:
% - 2005 Aug  7  - fpereira - added average pvalue under null hypothesis of no correlation with each neighbour
% - 2005 July 10 - fpereira - created
%
%

function [correlationImages,pvalueImages,meanImages,stdvImages,conditions] = summarizeIDM_voxelNeighbourCorrelation( varargin )

%
% Process parameters
%

this = 'summarizeIDM_voxelNeighbourCorrelation';

l = length(varargin);
if l < 3
  eval(sprintf('help %s;',this)); return;
else
  info = varargin{1};
  data = varargin{2};
  meta = varargin{3};
end

%% determine other information from them

trialConds = [info(:).cond]; nTrials = length(trialConds);
trialLens  = [info(:).len];
useTrials  = find(trialConds > 0); nUse = length(useTrials);
useLens    = trialLens(useTrials);
conditions = unique(trialConds(useTrials)); nConditions = length(conditions);
nVoxels    = size(data{1},2);

correlationImages = zeros(nConditions,nVoxels);
pvalueImages      = zeros(nConditions,nVoxels);
meanImages        = zeros(nConditions,nVoxels);
stdvImages        = zeros(nConditions,nVoxels);

% find columns of voxels that neighbour each voxel
fprintf('%s: finding neighbours of each voxel\n',this);
[columnsToNeighbours,numberOfNeighbours] = findColumnNeighbours(meta);

%
% Compute image for each condition;
%

for c = 1:nConditions
  
  cond = conditions(c);
  
  trials = find(trialConds == cond);
  lens   = trialLens(trials);
  total  = sum(lens);
  
  fprintf('%s: computing correlation image for condition %d from %d trials\n',this,cond,length(trials));
  
  %  Concatenate trials of interest, per condition
  matrix = zeros(total,nVoxels);
  
  idx = 1;
  for nt = 1:length(trials)
    
    trial = trials(nt);
    len   = lens(nt);
    
    matrix(idx:(idx+len-1),:) = data{trial};
    idx = idx + len;
  end

  %% compute correlation
  [correlationImages(c,:),pvalueImages(c,:),meanImage(c,:),stdvImage(c,:)] = compute_correlationWithNeighbours(matrix,columnsToNeighbours,numberOfNeighbours); 

end

%
% Compute correlation across time of a voxel with each of its neighbours
%
% In:
% - a #time points x #voxels matrix (1 trial in a Data structure of
% an IDM, say)
% - the Meta from the IDM
%
% Out:
% - a #voxels image of correlations
%
% History:
% 28 Apr 04 - fp - created

function [correlationImage,pvalueImage,meanImage,stdvImage] = compute_correlationWithNeighbours( matrix,columnsToNeighbours,numberOfNeighbours )

% compute statistics that will be needed
meanImage = mean(matrix,1);
stdvImage = std(matrix,1,1);
[len,nvoxels] = size(matrix);

% subtract the mean from the matrix
matrix = matrix - ones(len,1)*meanImage;

%% Finally compute average correlation with neighbours
correlationImage = zeros(nvoxels,1);
pvalueImage      = zeros(nvoxels,1); % not implemented yet

for v = 1:nvoxels
  
  % find neighbours of this voxel
  nn = numberOfNeighbours(v);
  if nn
    neighbours = columnsToNeighbours(v,1:nn);
  
    % compute correlation of voxel with all neighbours
    covars = repmat(matrix(:,v),1,nn) .* matrix(:,neighbours);
    covars = sum(covars,1);
    sigmas = stdvImage(neighbours) * stdvImage(v) * len;

    correlations = covars ./ sigmas;
    correlationImage(v) = mean(correlations);
    
    % Average p-value of all pairwise correlations.
    % The correlation between each two voxels can be given
    % a p-value under the null hypothesis that there is no
    % correlation.
    pvalues = compute_correlationScorePvalues( correlations, len );
    pvalueImage(v) = median(pvalues);
  else
    correlationImage(v) = 0;
  end

  % TODO: this shouldn't happen now, but while I'm debugging...
  if isnan(correlationImage(v))
    nn
    stdvImage(neighbours)
    correlations
    covars
    sigmas
    pause
  end
  
end

%
% Takes in a vector of correlation computation scores, each one being
% correlation between the time series of two voxels (length n). Computes the
% probability of each score under the null hypothesis that the two 
% voxels used are not correlated
%

function [pvalues] = compute_correlationScorePvalues( scores, n )

if n <= 5
  fprintf('warning: you are computing correlations between two very short time series (%d points)\n',n);
  pause; return;
end
  
normalValues = sqrt(n-3)/2 * ( log( (1+scores)./(1-scores) ) );
pvalues      = 1 - normcdf(abs(normalValues));


% Find the columns that contain neighbours of each voxel.
% Outputs:
% - columnsToNeighbours: #voxels x 26 (only the first #ofneighbours positions are filled)
% - # of neighbours: #voxels x 1 (number of neighbours for voxel in columnsToNeighbours)
%

function [columnsToNeighbours,numberOfNeighbours] = findColumnNeighbours(meta)    
  nvoxels = size(meta.colToCoord,1);
  columnsToNeighbours = zeros(nvoxels,26);
  numberOfNeighbours  = zeros(nvoxels,1);

  for v = 1:nvoxels
    tmp1 = meta.colToCoord; tmp1(v,:) = [1000 1000 1000];
    tmp2 = ones(nvoxels,1) * meta.colToCoord(v,:);
    tmp1 = (tmp1 - tmp2).^2;
    tmp1 = sum(tmp1,2); % neighbours have scores [0,3]
    neighbours = find(tmp1 <= 3);
    numberOfNeighbours(v) = length(neighbours);
    columnsToNeighbours(v,1:numberOfNeighbours(v)) = neighbours';
  end
  

function [] = quickPlot( imageToPlot, meta )

% Find the linearized index for each IDM column in a 3D matrix
dimx = meta.dimx; dimy = meta.dimy; dimz = meta.dimz; dimensions = [dimx dimy dimz];
tmp  = meta.colToCoord;
indicesIn3D = sub2ind(dimensions,tmp(:,1),tmp(:,2),tmp(:,3));

% each image will be placed in this volume
volume = repmat(NaN,[dimx,dimy,dimz]);
nrows = ceil(sqrt(dimz)); ncols = nrows;

%
% Plot the correlation images
%

scale = [-1 1];

volume(indicesIn3D) = imageToPlot;

for z = 1:dimz
  
  subplot(nrows,ncols,z);
  
  imagesc(rot90(volume(:,:,z),-1),scale);
  if z == 1; colorbar('vertical'); end
  view(0,90);
  
  set(gca,'XTick',[]);set(gca,'YTick',[]);axis square;
end


