% Computes t-tests of difference in mean signal between baseline and activity conditions
%
% In:
% - info+data+meta
% - (optional) test - an array where each row has two elements
% indicating a pair of conditions to test against each other. If
% absent, the program defaults to testing 1 vs all other conditions
%
% e.g. to test 1 vs 2 and 2 vs 3, use as argument the array
%      [1 2
%       2 3]
%
% The test performed for the pair [cond1,cond2] is
%
% mean(1) < mean(condition) if one of the conditions is 1
%
% or
%
% mean(cond1) != mean(cond2)
%
% if none is condition 1
%
% Out:
% - results - an upper triangular cell array with results of all
% the requested tests. I.e.

% test vs cond1 cond2 cond3 ... cond n
% cond 1        +     +         +  
% cond 2              +         +
% ...
% cond n-1                      +
% cond n
%
% where each cell contains a map of test p-values at each voxel.
%
% This is a vector where position v contains
%
% probability that the mean signal at this voxel is not different for the
% two conditions being tested
%
% Therefore, the smaller the value the more likely it is that the
% means are different.
%
% Dependencies:
%
% History
% - 020731 - fpereira - wrote and debugged with test_this below
% - 28 Oct 02 - fp - added comments and rechecked
%

function [results] = mri_computeTvalues ( varargin )
  
  % process arguments
  
  l = length(varargin);
  if l < 3
    fprintf('syntax: mri_computeTvalues(info,data,meta,[tests]])\n');
    fprintf(['- tests - array where each row has a pair of conditions' ...
	     ' to test against each other. Defaults to 1 vs all' ...
	      ' others. 0 denotes condition 1 vs non-1(combined).\n']);    
    return;
  else
    info = varargin{1};
    data = varargin{2};
    meta = varargin{3};
  end
  
  % WARNING: HACK to use the precomputed per-ROI t-maps instead of
  % building them again.
  % Works only for the brainlex study, for now

  study   = meta.study;
  subject = meta.subject;
  
  switch study
   case {'data-brainlex','data-syntamb2','data-categories'}
    
    if strcmpi(deblank(meta.roi),'CUSTOM')
      % just go ahead and compute everything
    else      
      % unpack the roi names
      R = deblank(meta.roi);
      nrois = 0;
      [T,R] = strtok(R,'_');      
      while ~isempty(T)
	nrois = nrois + 1;
	rois{nrois} = T;
	[T,R] = strtok(R,'_');
      end
  
      fprintf('Using precomputed activity t-test maps for ROIs\n');
      rois
      
      [results] = loadAndMergeROItmaps(subject,study,rois);
      return;
    end
    
   otherwise
    % keep going and compute everything
  end
    
  % figure out data details and run sanity checks
  %[ntrials,nvoxels,nconds,minTrialLenCond,ntrialsCond] = IDMinformation( info,data,meta);
  IDM_information=IDMinformation( info,data,meta,meta.study );
  ntrials=IDM_information.nTrials;
  nvoxels=IDM_information.nVoxels;
  nconds=IDM_information.nConds;
  minTrialLenCond=IDM_information.minTrialLenCond;
  ntrialsCond=IDM_information.nTrialsCond;
  
  
  if l > 3
    condPairsToTest = varargin{4};
  else
    condPairsToTest = [ones(nconds-1,1),(2:nconds)'];
  end

  % compile all the data into #conditions separate bins  
  databins = cell(nconds,1);
  binidx   = zeros(nconds,1);
  clen     = zeros(nconds,1); % total trial length cond c
  
  for c=1:1:nconds

    databins{c} = cell(ntrialsCond(c),1);

    for i=1:1:ntrials
      cond = info(i).cond;
      len  = info(i).len;
      if cond > 0 & (cond == c)
	binidx(c) = binidx(c) + 1;

	databins{c}{binidx(c)} = data{i};
	clen(c)   = clen(c) + len;
      end
    end
  end

  % merge each bin into a megatrial
  cdata = cell(nconds,1);

  for c=1:1:nconds
    cdata{c} = zeros(clen(c),nvoxels);
    idx = 1;
    
    for bi=1:1:binidx(c)
      nextidx = idx + size(databins{c}{bi},1);
      cdata{c}(idx:1:(nextidx-1),:) = databins{c}{bi}(:,:);
      idx = nextidx;
    end
  end

  % Run t-tests of the difference between the means of each pair of
  % conditions requested in the arguments.
  %
  % The t-test performed depends on whether one of [cond1,cond2] is
  % equal to 1. More specifically:
  %
  % mean(1) < mean(condition) if one of the conditions is 1
  %
  % or
  %
  % mean(cond1) != mean(cond2)
  %
  % if none is condition 1.
  %
  % The test done is given by the tail argument passed to ttest2
  
  ntests  = size(condPairsToTest,1);
  results = cell(nconds,nconds);
  alpha   = 0.01; % not used
  
  for t=1:ntests
    c1 = condPairsToTest(t,1);
    c2 = condPairsToTest(t,2);   

    fprintf('\tt-testing %d vs %d\n',c1,c2);
    
    % pick type of test (one-tailed, two-tailed)
    if c1==1; tail=-1; elseif c2==1; tail=1; else tail=0; end
    
    % I'm afraid I couldn't think of any way around doing this
    % voxel by voxel. The ttest2 function requires two vectors as
    % input
    pvalues = zeros(1,nvoxels);
    for v=1:nvoxels 
      [h,pvalue,ci] = ttest2( cdata{c1}(:,v), cdata{c2}(:,v), alpha, tail);  
      % something odd is happening, it seems it's returning 1-pvalue
      %      [h,pvalue] = kstest2( cdata{c1}(:,v), cdata{c2}(:,v), alpha, -1);  
      pvalues(v) = pvalue;
    end   
  
    results{c1,c2} = pvalues;
  end

  

% WARNING: function used to load precomputed results for this
% function
function [results] = hack_readPrecomputedResultsAll(subject,study)
  
  information = mri_reference(study);
  dataplace   = information.dataplace;
  location    = sprintf('%s/%s/activityTtests',dataplace,study);
  % fetch data
  load(sprintf('%s/%s.mat',location,subject));
  
function [results] = hack_readPrecomputedResultsROI(subject,study,roi)
  
  information = mri_reference(study);
  dataplace   = information.dataplace;
  location    = sprintf('%s/%s/activityTtests',dataplace,study);
  % fetch data
  load(sprintf('%s/ttest.%s.%s.%s.mat',location,subject,roi,study));;

  
function [results] = loadAndMergeROItmaps(subject,study,rois)
  
  information = localInformation(study);
  dataplace   = information.dataplace;
  location    = sprintf('%s/%s/%s/data',dataplace,study,subject);
  nrois       = length(rois);

  fprintf('\tloading precomputed tmaps\n');
  
  % load result maps + meta information for its ROI
  for r=1:nrois
    roi = rois{r};
    load(sprintf('%s/detrend.CCBIttest.%s.mat',location,roi));
    roiResults{r} = results;
    load(sprintf('%s/detrend.%s.mat',location,roi));
    meta.roi = roi;
    roiMetas{r}   = meta;
    nvoxels(r)    = size(meta.colToCoord,1);
  end
  clear info data meta results;
    
  nconds = size(roiResults{1},2);
  
  % pack the results into IDMs for easier merging
  % - one IDM per ROI
  % - each trial in the IDM contains the map for one condition
  for r=1:nrois
    data{r} = cell(nconds,1);
    for c=1:nconds    
      info{r}(c).cond = c;
      info{r}(c).len  = 1;
      results = roiResults{r};
      data{r}{c} = results{1,c};
    end
    data{r}{1} = zeros(1,nvoxels(r)); % need to have values
    meta{r} = roiMetas{r};
  end
    
  % now use an IDM merging tool to merge the t-map;
  argString = 'info{1},data{1},meta{1}';
  for r=2:nrois; argString = [argString,sprintf(',info{%d},data{%d},meta{%d}',r,r,r)]; end
  cmd = sprintf('[minfo,mdata,mmeta] = transformIDM_mergeMulti(%s);',argString);
  eval(cmd);
  
  % and reassemble a t-map result structure
  results = cell(nconds,nconds);
  for c = 2:nconds
    results{1,c} = mdata{c};
  end
  

function [] = test_this()
  
% do a mock ROI with 4 voxels, and 3 conditions.
% voxel 1 - inactive during all conditions
% voxel 2 - active during condition 2
% voxel 3 - active during condition 3
% voxel 4 - active during conditions 2 and 3
  
  pointsPerTrial = 50;
  meta.nvoxels   = 4;
  
  % 2 trials in each of conditions 1,2 and 3  
  tinfo(1).cond = 1;
  tinfo(1).len  = pointsPerTrial;
  tinfo(2).cond = 2;
  tinfo(2).len  = pointsPerTrial;
  tinfo(3).cond = 3;
  tinfo(3).len  = pointsPerTrial;
  tinfo(4).cond = 2;
  tinfo(4).len  = pointsPerTrial;
  tinfo(5).cond = 1;
  tinfo(5).len  = pointsPerTrial;
  tinfo(6).cond = 3;
  tinfo(6).len  = pointsPerTrial;
  tmeta.ntrials = 6;
  tmeta.study   = 'dummy';
  
  % means for each of four voxels in the three conditions
  % 
  means{1} = [0 0 0 0];
  stdvs{1} = [0.5 0.5 0.5 0.5];
  means{2} = [0 0 4 6];
  stdvs{2} = [0.5 0.5 0.5 0.5];
  means{3} = [0 4 0 6];
  stdvs{3} = [0.5 0.5 0.5 0.5];
  tdata = cell(tmeta.ntrials,1);
  
  % generate dataset
  for nt=1:tmeta.ntrials
    len  = tinfo(nt).len;
    cond = tinfo(nt).cond;
   
    tdata{nt} = randn(len,meta.nvoxels) ./ (ones(len,1)*stdvs{cond}) + (ones(len,1)*means{cond});    
  end
  tmaps = mri_computeTvalues(tinfo,tdata,tmeta);
