% plotIDM_sliceTimeSeries(info,data,meta,trial,slice)
%
% Plots a given slice of a brain volume. Inside each voxel
% is a plot of the timeseries of that voxel during a particular
% trial.
%
% Suggested use:
% - apply transformIDM_meanTrialCond.m to the IDM beforehand to
%   get an IDM with one average trial per condition
% - call this code on the average IDM and request the trial for the
%   condition you are interested in.
%
% Inputs:
% - IDM
% - trial to plot
% - slice to plot
%
% Outputs:
% - none, it just opens a figure with the desired plot
%
% History:
% - 2004 Oct 26 - fpereira - created
%
% Notes:
%

function [] = plotIDM_sliceTimeseries( varargin )

  DEBUG = 1;

  %% process arguments
  l = length(varargin);
  if l < 5; help plotIDM_sliceTimeseries; return;
  else
    info  = varargin{1};
    data  = varargin{2};
    meta  = varargin{3};
    trial = varargin{4};
    slice = varargin{5};
  end
    
  %% figure out things  
  nTrials     = length(data);
  trialConds  = [info(:).cond];
  trialLens   = [info(:).len];
  totalLength = sum(trialLens); 
  nVoxels     = size(data{1},2);
  fixTrials   = find(trialConds == 1);
  fixLens     = trialLens(fixTrials);
  dimx = meta.dimx; dimy = meta.dimy; dimz = meta.dimz; 
  
  % find coordinates of the voxels in this slice
  tmp     = meta.colToCoord;
  indices = find(tmp(:,3) == slice);
  coordinates3D = tmp(indices,:);
  coordinates2D = coordinates3D(:,[1 2]);
  indices2D = sub2ind([dimx dimy],coordinates2D(:,1),coordinates2D(:,2));
  indices3D = sub2ind([dimx dimy dimz],coordinates3D(:,1),coordinates3D(:,2),coordinates3D(:,3));
  coordToCol2D = zeros(dimx,dimy);
  coordToCol2D(indices2D) = 1:length(indices2D);
  
  % 2D mask - 1 if voxel(x,y) is in IDM, 0 otherwise
  mask2D = zeros(dimx,dimy);
  mask2D(indices2D) = 1;
  
  %
  % Plot
  %
  % loop over the 2 dimensions placing axes in the grid position if
  % the voxel is in the IDM
  
  % padding between panes, where units are "fraction of the size of a pane"
  vpaddingInterPair = 0.5;
  hpaddingInterPair = 0.5;
  
  % computations of distances, taking into account dimx,dimy and padding
  nvertical   = dimy + (dimy-1)*vpaddingInterPair;
  nhorizontal = dimx + floor((dimx-1)*hpaddingInterPair) + 1;
  spaceWidth  = 1/nhorizontal;
  spaceHeight = 1/nvertical;
  
  % computations pertaining to the # of axes within a pane
  axisWidth  = spaceWidth;
  axisWidth  = 0.99*axisWidth;
  axisHeight = spaceHeight;
  axisHeight = 0.99*axisHeight;

  hidx = 0;
  ncols = dimx;
  nrows = dimy;

  for c = 1:ncols
    
    vidx = 1 - spaceHeight;
    for r = 1:nrows

      x = c;
      y = r;
      
      if mask2D(x,y)
	% plot a time series in the axis for this voxel
	axes('position',[hidx vidx axisWidth axisHeight]);
	voxel = coordToCol2D(x,y);
	plot( data{trial}(:,voxel) );
	set(gca,'XTick',[]); set(gca,'YTick',[]);
      end
	
      vidx = vidx - spaceHeight - vpaddingInterPair*spaceHeight;
    end
    
    hidx = hidx + spaceWidth + hpaddingInterPair*spaceWidth;  
  end
  
  return
