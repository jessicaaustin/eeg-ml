% Create a new IDM containing only the selected subset of voxels
%
% In:
% - info, data and meta
% - a list of voxels
%
% Out:
% - info, data, meta restricted to selected voxels in order given
%
% Dep:
%
% Notes
% - 14 Mar 2005 - fpereira - adapted from previous versions

function [info,data,meta] = transformIDM_selectVoxelSubset( varargin )

  l = length(varargin);

  if l < 4
    fprintf('syntax: transformIDM_selectVoxelSubset(i,d,m,list of voxel #s)\n');
    return;
  else  
    info = varargin{1};
    data = varargin{2};
    meta = varargin{3};
    voxels = varargin{4};
  end

  ntrials = length(data);
  
  for i=1:1:ntrials    
    data{i} = data{i}(:,voxels);
  end
    
  meta.nvoxels    = length(voxels);
  meta.colToCoord = meta.colToCoord(voxels,:);
  meta.coordToCol = zeros(meta.dimx,meta.dimy,meta.dimz);
  xcoords = meta.colToCoord(:,1);
  ycoords = meta.colToCoord(:,2);
  zcoords = meta.colToCoord(:,3);
  indices = sub2ind([meta.dimx meta.dimy meta.dimz],xcoords,ycoords,zcoords);  
  meta.coordToCol(indices) = 1:meta.nvoxels;
  
  
function test_this()
  
  a = magic(7);
  
  actualmean = mean(a,2);
  actualstd  = std(a,0,2);
  
  ginfo(1).len = 7;
  ginfo(1).cond = 1;
  gmeta = 1;
  gdata{1} = a;
  gmeta.study = 'data-categories';  
  gmeta.dimx=64; gmeta.dimy=64; gmeta.dimz=16;
  gmeta.coordToCol = zeros(64,64,16);

  for v=1:7
    gmeta.colToCoord(v,:) = [1 1 v];
    gmeta.coordToCol(1,1,v) = v;
  end
    
  [ia,da,sa] = transformIDM_selectVoxelSubset(ginfo,gdata,gmeta);
