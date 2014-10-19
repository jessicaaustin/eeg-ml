% Plot the mean barin image according to each slice
%
% Inputs:
% - info, data, meta
%
% History:
% - 2005 May 06 - Wei 
% adapted from FP's code
%
% Notes:
%

function plotBrain( varargin )

  %% process arguments
  l = length(varargin);
  if l < 3; help plotBrain; return;
  else
    info = varargin{1};
    data = varargin{2};
    meta = varargin{3};
  end
  if l ==4
      threshold = varargin{4};
  else
      threshold = -100000000;
  end
  
  %% compute a mean image - avgDataMatrix
  
  nTrials     = length(data);
  trialConds  = [info(:).cond];
  trialLens   = [info(:).len];
  totalLength = sum(trialLens); 
  nVoxels     = size(data{1},2);

  avgDataMatrix = zeros(1,nVoxels);    
  for nt = 1:nTrials; avgDataMatrix = avgDataMatrix + sum(data{nt},1); end
  avgDataMatrix = avgDataMatrix / totalLength;
  
  %% threshold and update meta and data

  % find voxels above threshold;
  mask = avgDataMatrix > threshold;
  voxelsToKeep  = find(mask);
  nvoxelsToKeep = length(voxelsToKeep);
  
  % update meta;
  dimx = meta.dimx; dimy = meta.dimy; dimz = meta.dimz;
  meta.colToCoord = meta.colToCoord(voxelsToKeep,:);
  tmp = meta.colToCoord;
  newIndicesIn3D = sub2ind([dimx dimy dimz],tmp(:,1),tmp(:,2),tmp(:,3));  
  meta.coordToCol = zeros(dimx,dimy,dimz);
  meta.coordToCol(newIndicesIn3D) = 1:nvoxelsToKeep;
  
  % update data;
  for nt = 1:nTrials    
    data{nt} = data{nt}(:,voxelsToKeep);
  end

  


  % find voxel positions in 3D
  tmp  = meta.colToCoord;
  indicesIn3D = sub2ind([dimx dimy dimz],tmp(:,1),tmp(:,2),tmp(:,3));
    
  % place their average values in a 3D volume
  volume = repmat(NaN,[dimx,dimy,dimz]); % buffer to keep volumes read
  volume(indicesIn3D) = avgDataMatrix(voxelsToKeep);
    
  % plot
  scale = [prctile(avgDataMatrix(voxelsToKeep),1) prctile(avgDataMatrix(voxelsToKeep),99)];
  nrows = ceil(sqrt(dimz)); ncols = nrows;
  
  for z = 1:dimz	  
    subplot(nrows,ncols,z);      
    imagesc(volume(:,:,z),scale);
  end    


