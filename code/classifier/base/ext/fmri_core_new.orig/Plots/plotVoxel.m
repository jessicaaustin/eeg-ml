%  Plot the time course for the selected voxel, concatenating all
%  trials.  Superimpose a plot of condition numbers.
%
% plotVoxel(info,data,meta,x,y,z,<titlestring>)
%   x,y,z : coordinates of the voxel to be plotted
%   titlestring : optional text string to be printed with plots.
%
% or
%
% plotVoxel(info,data,meta,<column number>,<titlestring>)
%
% Example:
%   plotVoxel(i,d,m,46,54,10)  % plots time course for voxel <46,54,10>
%  
% History: 
% 8/16/02 tom - created. (based on plotTrial)
% 10/16/02 fp - adapted to accept a column number as an argument

function [] = plotVoxel( varargin )

  % process arguments
  l = length(varargin);
  if l < 4
    fprintf('syntax:\nplotVoxel(<info>,<data>,<meta>,<x>,<y>,<z>,[titlestring])\nor\nplotVoxel(<info>,<data>,<meta>,<voxel column>,[titlestring])\n'); 
    return;
  end

  info  = varargin{1};
  data  = varargin{2};
  meta  = varargin{3};
  titlestring   = '';
  
  if l > 5  
    x = varargin{4};
    y = varargin{5};
    z = varargin{6};  
    if l > 6
      titlestring = varargin{7};
    end
  else
    column = varargin{4};
    if l > 4
      titlestring = varargin{5};
    end
    
    x = meta.colToCoord(column,1);
    y = meta.colToCoord(column,2);
    z = meta.colToCoord(column,3);
  end
    
  % delete graphics in the current figure window
  clf reset;
  
  vtc = getVoxelTimecourse(info,data,meta,x,y,z);
  ctc = getConditionTimecourse(info,data,meta);

  tstr=sprintf('%s %s subject%s, region %s,', titlestring,meta.study,meta.subject,meta.roi);
  tstr=sprintf('%s\n <x,y,z>=<%d,%d,%d> ', tstr,x,y,z);
  conditionstr= sprintf(' %0.5g', [info.cond]);
  tstr=sprintf('%s condition sequence = %s', tstr, conditionstr);
  
  plot(vtc);
  hold on;
  plot(ctc, 'r');
  title(tstr);
  
  

