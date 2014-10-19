% Merge two or more IDMs (allows overlapping voxels)
%
% In:
% - info, data and meta for each
%
% Out:
% - the joint info, data and meta
%
% Notes:
% - The token information must match (#trials, lengths and conditions)
%

function [minfo,mdata,mmeta] =transformIDM_mergeMulti( varargin );
  
  l = length(varargin);
  
  if mod(l,3) | l < 3
    fprintf('syntax: transformIDM_mergeMulti(info1,data1,meta1,info2,data2,meta2,[infok,datak,metak])\n');
    return;
  else
    nIDMs  = l/3;
    cursor = 1;
    
    for i=1:nIDMs
      info{i} = varargin{cursor}; cursor = cursor + 1;
      data{i} = varargin{cursor}; cursor = cursor + 1;
      meta{i} = varargin{cursor}; cursor = cursor + 1;      
    end
    
    if nIDMs == 3
      % no need to merge
      minfo = info; mdata = data; mmeta = meta; return;
    end  
  end
  
  % sanity checks
  % # trials, their lengths and conditions must match
  
  % assuming they all do...
  
  % Start creating the merged IDM
  ntrials      = length(data{1});
  mmeta.ntrials = ntrials;
  minfo = info{1};
  mdata = cell(ntrials,1);
  mmeta.dimx = meta{1}.dimx;
  mmeta.dimy = meta{1}.dimy;
  mmeta.dimz = meta{1}.dimz;
  mergedROI = [];
  
  % merge data

  % 1) get indices of voxels in a vectorization of the 3D matrix
  totalNvoxels = 0;
  for i=1:nIDMs  
    ix = meta{i}.colToCoord(:,1);
    iy = meta{i}.colToCoord(:,2);
    iz = meta{i}.colToCoord(:,3);
    indices{i} = (sub2ind([mmeta.dimx,mmeta.dimy,mmeta.dimz],ix,iy,iz))';
    nvoxels(i) = length(indices{i});
    totalNvoxels = totalNvoxels + nvoxels(i);
    if isempty(mergedROI)
      mergedROI = meta{i}.roi;
    else;
      mergedROI = sprintf('%s_%s',mergedROI,meta{i}.roi);
    end
  end
  
  % 2) merge them and convert back to subscripts over the 3D matrix
  allIndices = zeros(1,totalNvoxels);
  cursor     = 1;
  for i=1:nIDMs
    allIndices(cursor:(cursor+nvoxels(i)-1)) = indices{i};
    cursor = cursor + nvoxels(i);
  end

  allIndices   = unique(allIndices);
  totalNvoxels = length(allIndices); % there may be overlap 
  [subx,suby,subz] = ind2sub([mmeta.dimx,mmeta.dimy,mmeta.dimz],allIndices);

  % 3) finally merge the data into matrix m, trial by trial
  mmeta.nvoxels = totalNvoxels;
  
  m = zeros(mmeta.dimx,mmeta.dimy,mmeta.dimz);
    
  for nt=1:1:ntrials
    for t=1:minfo(nt).len;
%      fprintf('\ttrial %d time %d\n',nt,t);
      for i=1:nIDMs
	% place data in the appropriate location
	% of a common matrix
	m(indices{i}) = data{i}{nt}(t,:);
      end
      
      % dump the common matrix into a trial      
      mdata{nt}(t,:) = m(allIndices);
    end
  end
  
  mmeta.colToCoord = [subx',suby',subz'];  
  mmeta.coordToCol = zeros(mmeta.dimx,mmeta.dimy,mmeta.dimz);
  mmeta.coordToCol(allIndices) = 1:mmeta.nvoxels;
  mmeta.roi = mergedROI;
    

  
  
function [] = testThis()
  
  a = magic(4);
  
  cube = zeros(4,4,3);
  cube(:,:,1) = a;
  cube(:,:,2) = a + 20;
  cube(:,:,3) = a + 40;

  trial = zeros(10,4*4*3);
  
  for t=1:1:10
    for z=1:1:3
      for y=1:1:4
	for x=1:1:4
	  col = x + (y-1)*4 + (z-1)*4*4;
	  trial(t,col) = cube(x,y,z)+1;
	end
      end
    end
  end
  

  
  % split it into four parts
