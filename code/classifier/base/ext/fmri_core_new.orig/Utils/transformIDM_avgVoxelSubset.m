% Average across a group of voxels
%
% In:
% - info, data and meta
% - optionally, a list of columns [default all] 
% - optionally, a 1*n weight matrix for the column list [default equal weight] 

% Out:
% - an IDM with a single voxel, the average of the subset given
% - standard deviation vector (parallels data through time)
%
% Dep:
%
% Notes:
% 01 Nov 02 - fp -created

function [info,avgdata,meta,stdev] = transformIDM_avgVoxelSubset ( varargin )

  l = length(varargin);

  if l < 3
    fprintf('syntax: transformIDM_avgVoxelSubset(i,d,m,[column list,[weight list]]\n');
    return;
  else  
    info = varargin{1};
    data = varargin{2};
    meta = varargin{3};
    meta.nvoxels = 1; % average will have only one
    
    ntrials = max(size(info));
    nvoxels = size(data{1},2);
    
    if l > 3  
      columnlist = varargin{4};
      if l > 4
	mixing    = varargin{5};  
	normalAvg = 0;
      else
	mixing = ones(size(columnlist));
	normalAvg = 1;
      end
    else
      columnlist = 1:1:nvoxels;
      mixing     = ones(1,nvoxels);
      normalAvg  = 1;
    end
    
    if ~isequal(size(mixing),size(columnlist))
      fprintf(1,'error: transformIDM_avgVoxelSubject: the size of the mixing matrix should match that of the columnlist');
      return;
    end
  end
    
  mixing = mixing ./ sum(mixing); % make mixing weights add up to 1
      
  avgdata = cell(ntrials,1);
  stdev   = cell(ntrials,1);
  tmp     = cell(ntrials,1);
  nvoxels = size(data{1},2);
  ncols   = length(columnlist);
  
  for i=1:1:ntrials
    len = info(i).len;
    
    if l <= 3
      % average of all cols, without weight matrix
      % spare copying of everything
      avgdata{i} = mean(data{i},2);
      stdev{i}   = std(data{i},0,2);
    else
      if normalAvg
	avgdata{i} = mean(data{i}(:,columnlist),2);
	stdev{i}   = std( data{i}(:,columnlist),0,2);
      else
	% mean applying mixing weights
	tmp = zeros(len,ncols);
	tmp(:,:) = data{i}(:,columnlist); 
	avgdata{i} = sum(tmp .* (ones(1,len)' * mixing),2);
	tmp = (tmp - avgdata{i} * ones(1,ncols)).^2; % (x_i - mean)^2
	tmp = tmp .* (ones(1,len)' * mixing); % weighed by mixing matrix
	stdev{i}   = sqrt(sum( tmp, 2));
      end
    end
  end
    
function test_this()
  
  a = magic(7);
  
  actualmean = mean(a,2);
  actualstd  = std(a,0,2);
  
  ginfo(1).len = 7;
  gmeta = 1;
  gdata{1} = a;
  
  [ia,da,sa] = transformIDM_avgVoxelSubset(ginfo,gdata,gmeta);
  
  fprintf(1,'should all be 1 m=%d, s=%d\n',isequal(da{1},actualmean),isequal(sa{1},actualstd));