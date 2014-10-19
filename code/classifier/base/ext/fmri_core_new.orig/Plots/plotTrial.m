%  Plot the voxel time courses for trial number trialNum, one z slice at a
%  time, for the given IDM. Lays out voxel plots according to their geometric
%  relationship.
%
% plotTrial(info,data,meta,trialNum,[titlestring],[dointeractive=0|1 (default=1)])
%
% Voxel plots are laid out in same spatial relationship as voxels
% themselves, one z slice at a time.
% trialNum : number of the trial from data that is to be displayed
% titlestring : an optional text string to be printed with plots.
% dointeractive : if 1, display to terminal, else print to file
%  
% History: 
% 1/31/02 Tom - created.  Next, add auto scaling for min,max Activation
% 03/06/02: fp - modified this to make it a module callable by a
% wrapper that iterates over all subjects and ROIs of interest, and
% also tweaked it a little.
% 8/13/02 tom - remamed to plotTrial, cut out the hidden smoothing
% it did on signal before displaying it.  Now what you see is what
% the program sees.
%
% 4/11/05 tom - fixed a bug which occured if number of columns was odd
%

function [] = plotTrial( varargin )

  % process arguments
  l = length(varargin);
  if l < 4
    fprintf(1,'syntax: plotTrial(<info>,<data>,<meta>,<trial>,[titlestring],[dointeractive=0|1 (default=1)])\n'); 
    return;
  end

  info  = varargin{1};
  data  = varargin{2};
  meta  = varargin{3};
  trialNum = varargin{4};

  titlestring   = '';
  doInteractive = 1;
  if l > 4
    titlestring = varargin{5};
    if l > 5      
      doInteractive = varargin{6};
    end
  end
    
% delete current figure graphics 
  clf reset;

% parameters
  nvoxels      = size(data{1},2);
  colours = {'w'};
  assignments = ones(nvoxels,1);
    
% data dependent parameters
  trialBegin    = 1;
  trialEnd      = info(trialNum).len;
%  minActivation = min(min(data{trialNum}(trialBegin:1:trialEnd,:)));
%  maxActivation = max(max(data{trialNum}(trialBegin:1:trialEnd,:)));
  vmin = sort(min(data{trialNum}(trialBegin:1:trialEnd,:)));
  vmax = sort(max(data{trialNum}(trialBegin:1:trialEnd,:)));
  ntouse = floor(0.05*nvoxels);
  maxActivation = mean(vmax(nvoxels-ntouse+1:1:nvoxels));
  ntouse = floor(0.1*nvoxels);
  minActivation = mean(vmin(1:1:ntouse));
  fprintf(1,'min=%1.2f\tmax=%1.2f\n',minActivation,maxActivation);
%  return
  
  % collect the appropriate trial
  trialdata = data{trialNum};
  len       = info(trialNum).len;
  condition = info(trialNum).cond;
  
  % calculate maximum width and height of plot grid
  xMin    = min(meta.colToCoord(:,1));
  xMax    = max(meta.colToCoord(:,1));
  columns = 1 + xMax - xMin;
  yMin    = min(meta.colToCoord(:,2));
  yMax    = max(meta.colToCoord(:,2));
  rows    = 2 + yMax - yMin;   	% one extra to leave room for grid title
  xOffset = xMin-1;
  yOffset = yMin-1;
  disp([ xMin xMax yMin yMax rows columns]);

  assignments = ones(nvoxels,1);
  slices = unique(meta.colToCoord(:,3));
  nsubs  = rows * columns;
  
  for s=1:1:length(slices)
    z = slices(s);
    disp(z)
    
    % get row numbers in colToCoord for voxels in this slice
    sliceVoxels = find(meta.colToCoord(:,3)==z);

    % for each voxel in plane z, plot its time course at its xy
    % coordinate
    xminPlot=9999; xmaxPlot=-9999; yminPlot=9999; ymaxPlot=-9999;
    for v=1:1:length(sliceVoxels)
      coord = meta.colToCoord(sliceVoxels(v),:); 
      x=coord(1); y= coord(2);
      xminPlot=min(xminPlot,x); xmaxPlot=max(xmaxPlot,x);
      yminPlot=min(yminPlot,y); ymaxPlot=max(ymaxPlot,y);
      
	      
      timeCourse = trialdata(:,sliceVoxels(v));
      timeCourse = timeCourse( trialBegin:1:trialEnd,: );
            
      %      plotIdx = x-xOffset+((y-yOffset)*columns));
      plotIdx = nsubs - ((y-yOffset-1)*columns + (columns-(x-xOffset)));
      subplot(rows,columns,plotIdx);
      h = plot(timeCourse);	
      
      set(gca,'XTickLabel',{''});
      set(gca,'YTickLabel',{''});

      axis([0,trialEnd,minActivation+1,maxActivation]);      
    end;

        
    subplot(rows,columns,ceil(columns/2));
    set(gca,'XTickLabel',{''});
    set(gca,'YTickLabel',{''});
%    minActivation=-1;
%    maxActivation=5;
    
    tstr=sprintf('%s %s subject%s, trial number %d', titlestring,meta.study,meta.subject,trialNum);
    tstr=sprintf('%s\nregion %s\n',tstr,meta.roi);
    tstr=sprintf('%sz=%d x=[%d,%d] ',tstr,z,xminPlot,xmaxPlot);
    tstr=sprintf('%s y=[%d,%d], timeSteps=%d, condition=%d ', tstr,yminPlot,ymaxPlot,trialEnd,condition);
    tstr=sprintf('%s amplitude=[%1.1f %1.1f]', tstr, minActivation+1,maxActivation);
    title(tstr);
    
    if doInteractive
      fprintf(' press <return> to display next z slice');
      pause; clf;
    else
      outfile = sprintf('%ssubj%s_trial%d_z%d',titlestring,subject,trialNum,z);
      print('-dps',outfile);
    end
  end








