% loadSubjectdataMult(Study,Subject,ROIs,[separateROIs (default=1)])      
%
% Loads info+data+meta for specified study, subject, and ROIs from disk.
% Subject is ID of subject.  ROIs is a cell array of ROI names.  If the
% optional argument separateROIs? is 1, then create a distinct IDM for each
% ROI, otherwise merge all ROIs into a single IDM.
%
% Returns three cell arrays, one containing infos, one of meta's, one of
% data's, corresponding to the data loaded forr this study and subject
%
% Example:  
%  [is,ds,ms] = subjectLoadData('data-categories', '05502', {'LIT' 'CALC'}) 
%  sets is, ds, and ms to cell arrays, each containing two i's, d's, and
%  m's.  The first i,d, and m contain data for the ROI 'LIT', the second
%  contains data for the ROI 'CALC'.
%
% Dependencies:
% - process_dataToExamples depends on
%
% History
% - 020404 - fp
% - Added support for the categories data set and made example code
% standalone (rather than functions callable inside this)
% - Aug13,2002 Tom - added fields to meta: subject, study, roi
% - Aug18,2002 Tom - made separateROIs an optional argument
% - Sep 26,2002 Tom - edited filename generation to support windows or
%                                   as well as linux

function [rinfo,rdata,rmeta] = loadSubjectdataMult( varargin )
  present = pwd;
  l = length(varargin);
  
  if l < 1
    fprintf(1,'syntax: loadSubjectdata(study,subject,rois,[<use separate rois>],[preProcessSteps])\n\n');
    fprintf('- study - data-starplus-sp/ps or data-categories\n');
    fprintf('- subject - which subject inside the study\n');
    fprintf('- rois - which ROIs to use\n');
    fprintf('- useSeparateROIs (default=0) - if 1, return separate idm for each roi, else merge them all into a single idm\n');
    fprintf('- preProcessSteps - pre-processing to run on all the ROIs individually, regardless of whether they will be merged\n');
    fprintf('\n');
    return;
  else
    study           = varargin{1};
    subject         = varargin{2};
    rois            = varargin{3};
    useSeparateROIs = 1;
    doPreProcess    = 0;
    if l > 3
      useSeparateROIs = varargin{4};
      if l > 4
	preProcessSteps = varargin{5};
	doPreProcess    = 1;
      end	
    end
    nrois           = length(rois);
  end

  %
  % parameters
  ofd = 1;
  
  % collect relevant information
  information = localInformation(study);
  dataplace   = information.dataplace;
  host        = information.host;
  dimensions  = information.dimensions;
  dimx        = dimensions(1); dimy = dimensions(2); dimz = dimensions(3);
  nconds      = information.nconds;
  dataloc     = fullfile(dataplace,study);
    
  %
  % Start processing
  %
  
  if useSeparateROIs group='separate'; else group='joint'; end

  % 1) retrieve IDM for each ROI
  rinfo = cell(nrois,1);
  rdata = cell(nrois,1);
  rmeta = cell(nrois,1);
  dotherDataInfo = cell(nrois,1);
    
  for r=1:1:nrois

    roi = rois{r};
    info = []; % needed because info is the name of a matlab function
    roiDataFile = fullfile(dataloc,subject,'data',sprintf('detrend.%s.mat',roi));
    load(roiDataFile); % gives us info,data,meta
  
    [rinfo{r},rdata{r},rmeta{r}] = newFetchData( dataloc,subject,roi,study );
    rmeta{r}.roi   = roi;
    rmeta{r}.study = study;
    rmeta{r}.subject = subject;
    msg = sprintf('\t\tloading %s\t%d voxels\n',roi,rmeta{r}.nvoxels); fprintf(ofd,msg);
  
    if doPreProcess
      [rinfo{r},rdata{r},rmeta{r}] = transformIDM_preProcess(rinfo{r},rdata{r},rmeta{r},preProcessSteps);
    end
  
  end

  % 2) if ROIs are to be combined, merge them 
  % at the end we have a list of ROIs, useROIs, which either
  % contains a combined ROI or all of the ROIs
  
  if ~useSeparateROIs
    % gives us info,data,meta for first roi
    info = rinfo{1}; data = rdata{1}; meta = rmeta{1};
    msg = sprintf('\t\tmerging %s\n',rois{1}); fprintf(ofd,msg);

    % now merge left deep
    for r=1:1:nrois-1
      roi = rois{r+1};
      msg = sprintf('\t\tmerging %s\n',roi); fprintf(ofd,msg);
	
      % merge wih previous one
      [info data meta] = mri_mergeIDM( info,data,meta,rinfo{r+1},rdata{r+1},rmeta{r+1} );
    end	  
    clear rinfo rdata rmeta;
    
    rinfo{1} = info;
    rdata{1} = data;
    rmeta{1} = meta;
    clear info data meta;
    
    % build name for the joint ROI
    comboROI = rois{1};
    for r=2:1:nrois
      comboROI = sprintf('%s_%s',comboROI,rois{r});
    end
    useROIs{1} = comboROI;
    nuserois = 1;
    rmeta{1}.roi=comboROI;    
  else
    useROIs  = rois;
    nuserois = nrois;
  end
  
  for i=1:1:length(rmeta)
    rmeta{i}.study   = study;
    rmeta{i}.subject = subject;
    % CHECK with tom: This will assign the info and data into fields of the meta,
    % not pass them by reference
    %    rmeta{i}.info    = rinfo{1};
    %    rmeta{i}.data    = rdata{1};
    rmeta{i}.nsnapshots = sum([rinfo{i}.len]);
  end
  
  

%
% Capsule for the code that retrieves the data for a given data type/subject/roi
%

function [info,data,meta] = newFetchData( dataloc, subject, roi, study )

  info = []; % needed because info is also the name of a matlab function

  roiDataFile = fullfile(dataloc,subject,'data',sprintf('detrend.%s.mat',roi));
  load(roiDataFile); % gives us info,data,meta
  

% ATTENTION TOM: data type full will pluck out the simplest
% possible IDM - the rest is included so that you can get the code
% out for various IDM->IDM transformations, such as unfolding a
% block, averaging neighbouring (in time) images or picking active voxels

% - The code has to produce a valid info,data,meta.
% - Distilled data types are handled by the otherwise part of the
% switch, unless mentioned in the case statements
% - The case statements are meant to handle special formats such as
% edits of the data (voxel subsets, flattening of time, etc).

function [info,data,meta] = fetch_data( dataloc, subject, roi, dataType, study )
  
  % first get the data and do some preprocessing that is common to
  % a number of datatypes
  switch dataType
   
   case{'snapshot','full','abbrev','flat','unfolded','unfoldedWindowAvg'}
    info = []; % needed because info is the name of a matlab function
    roiDataFile = sprintf('%s/%s/data/detrend.%s.mat',dataloc,subject,roi);
    load(roiDataFile); % gives us info,data,meta
    
   case{'unfoldedSeparateBlocks'}
    % turn trials with the same condition into separate new conditions
    info = []; % needed because info is the name of a matlab function
    roiDataFile = sprintf('%s/%s/data/detrend.%s.mat',dataloc,subject,roi);
    load(roiDataFile); % gives us info,data,meta
    
    [info,data,meta] = process_separateByBlock(info,data,meta);
    
   case{'unfoldedAvgTrial'}
    % average all the trials in the same condition into a single one
    info = []; % needed because info is the name of a matlab function
    roiDataFile = sprintf('%s/%s/data/detrend.%s.mat',dataloc,subject,roi);
    load(roiDataFile); % gives us info,data,meta
    
    [info,data,meta,avgstd] = mri_averageIntraVoxel ( info, data, meta );

   case{'active','unfoldedActive','unfoldedActiveSeparateBlocks'}
    info = []; % needed because info is the name of a matlab function
    roiDataFile = sprintf('%s/%s/data/detrend.%s.mat',dataloc,subject,roi);
    load(roiDataFile); % gives us info,data,meta
    
    if isequal(dataType,'unfoldedActiveSeparateBlocks')
      [info,data,meta] = process_separateByBlock(info,data,meta);
    end
    
    % finds the top <nToAvg> most active voxels in each condition
    % and keeps all of those - at most #ROIS*#conditions*#nToKeep voxels
    % will be kept
    nToKeep = 10;
    nvoxels = size(data{1},2); % WARNING: fails if data{1}==[]
    activeCount = zeros(1,nvoxels);
    
    % get the Tvalues
    [results] = mri_computeTvalues (info,data,meta);
    nconds    = size(results,1);
    
    % calculate per condition ranks of voxels    
    for c=2:nconds
      pvalues   = results{1,c};
      ensemble  = sortrows(([pvalues;(1:nvoxels)])',1);
      topVoxels = ensemble(1:nToKeep,2);
      disp(topVoxels')
      activeCount(topVoxels) = activeCount(topVoxels) + 1;
%      fprintf('condition %d\n',c);
%      fprintf('\t%1.5f %d\n',(ensemble(1:nToKeep,:))');
    end

    activeVoxels = find(activeCount > 0); % these will be used
    nactive      = length(activeVoxels);
    
    % go through data and crop to these voxels
    % and update the meta information to reflect it
    nmeta = meta;
    nmeta.nvoxels = nactive;
    nmeta.colToCoord = meta.colToCoord(activeVoxels,:);
    nmeta.coordToCol = zeros(size(meta.coordToCol));
    % could be made more efficient, but few voxels anyway
    for v=1:nactive
      vc = nmeta.colToCoord(v,:);
      vx = vc(1); vy = vc(2); vz = vc(3);      
      nmeta.coordToCol(vx,vy,vz) = v;
    end
    
    ndata = cell(meta.ntrials,1);    
    for nt = 1:1:meta.ntrials
      ndata{nt} = data{nt}(:,activeVoxels);
    end
    
    data = ndata;
    meta = nmeta;
        
   otherwise
    % this is a distilled data abstraction
    % which has been pre computed
    roiDataFile = sprintf('%s/%s/data_distilled/distilled.%s.%s.mat',dataloc,subject,dataType,roi);
    load(roiDataFile);
  end
     
  % gather information about number,type and length of trials and
  % other study parameters
  [ntrials,nvoxels,nconds,minTrialLenCond,ntrialsCond] = mri_infoTrials( info,data,meta,study );
  
  information = mri_reference(study);
  dimensions  = information.dimensions;
  dimx        = dimensions(1); dimy = dimensions(2); dimz = dimensions(3);
  
  % decide how much of a trial to use - currently only used by
  % unfolded types

  switch study
   case {'data-starplus-sp','data-starplus-ps'}
    trialBegin = ones(nconds,1) + 6; % skip first 6 images
    trialEnd   = min(minTrialLenCond',ones(nconds,1) + 30);
   case {'data-categories'}
    trialBegin = ones(nconds,1) + 6; % skip first 3 images
    trialEnd   = min(minTrialLenCond',ones(nconds,1) + 37);  
   otherwise
    fprintf('WARNING: no support for unfolded syntamb2 studys\n');
  end
    
  % do any postprocessing  

  switch dataType
   case{'snapshot'}
    % full data set - all voxels, all trials, all time points in 1 trial
    [info,data,meta] = mri_snapshot(tSlice,info,data,meta);
    
   case{'full'}
    % full data set - all voxels, all trials, all time points in 1 trial
    
   case{'abbrev'}
    % finds the <nabbrev> voxels with higher avg signal value through time
    % and clips the data to comprise only that subset
    nabbrev = 20;

    % first identify the top x voxels in terms of activity
    [einfo,edata,emeta,es] = mri_averageIntraVoxel(info,data,meta);	    
    minlen = min(size(edata{2},1),size(edata{3},1));
    tdata  = edata{2}(1:1:minlen,:) + edata{3}(1:1:minlen,:);	    
    mtdata = mean(tdata,1);
    l      = length(mtdata);
    mtidx  = 1:1:l;
    sdata  = [mtidx',mtdata'];
    sdata  = sortrows(sdata,2);
    topcols = sdata(l-(nabbrev-1):1:l,1); % pick columns to use;
    topdata = sdata(l-(nabbrev-1):1:l,2); % pick columns to use;
				 %    [info,data,meta] =
				 %    mri_averageInterVoxel(info,data,meta,topcols');	 
    % then update the meta information to reflect it
    ndata = cell(meta.ntrials,1);
    nmeta = meta;
    for nt = 1:1:meta.ntrials
      ndata{nt} = data{nt}(:,topcols');
    end
    data = ndata;
    nmeta.nvoxels = nabbrev;
    nmeta.colToCoord = meta.colToCoord(topcols,:);
    nmeta.coordToCol = zeros(size(meta.coordToCol));

    for v=1:1:nabbrev
      x=nmeta.colToCoord(v,1); y=nmeta.colToCoord(v,2); z=nmeta.colToCoord(v,3);
      nmeta.coordToCol(x,y,z) = meta.coordToCol(x,y,z);
    end
    meta = nmeta;
        
   case{'active'}
      
   case{'flat'}
    % now flatten each trial to a single time point
    [info,data,meta] = mri_flattenTrialFilter(info,data,meta);

   
   case{'unfolded','unfoldedSeparateBlocks','unfoldedAvgTrial', ...
	'unfoldedActive','unfoldedActiveSeparateBlocks'}
    % takes each trial and turns each time point into a separate
    % trial of length 1 - this is useful for block designs
    % where the time doesn't matter (signal level should be the
    % same through most of the trial)
    
    % expand them    
    untrials = 1;
    
    % instrument to give us trial and example counts
    trialCounts = zeros(nconds,1);
    imgCounts   = zeros(nconds,1);

    for nt=1:1:ntrials
            
      len   = info(nt).len;
      cond  = info(nt).cond;
      tinfo = info(nt);
      tinfo.len = 1;

      if cond > 0
	s = trialBegin(cond); e = trialEnd(cond);
      	trialCounts(cond) = trialCounts(cond) + 1;
	imgCounts(cond)   = imgCounts(cond) + (e-s+1);
      else
	s = 1; e = len;
      end

      for t=s:1:e
	udata{untrials} = data{nt}(t,:);
	uinfo(untrials) = tinfo;
	untrials = untrials + 1;
      end
    end

    for c=1:nconds; fprintf('condition %d\t%d trials\t%d images\n',c,trialCounts(c),imgCounts(c));end

    data = udata;
    info = uinfo;
    meta.ntrials = untrials-1;

    % info+data+meta were transformed into a pseudo-dataset
    % with as many trials as images (minus the ones excluded)
    % and where each trial has length 1   

    
   case{'unfoldedWindowAvg'}
    % takes each trial and turns each time point into a separate
    % trial of length 1 - this is useful for block designs
    % where the time doesn't matter (signal level should be the
    % same through most of the trial)

    % Differs from "unfolded" in that each of these points is
    % replaced by the average of itself and a given number of
    % neighbours in both sizes(windowSize=1+#neighbours total);
    windowSize = 3;
    radius     = floor(windowSize/2);
            
    % expand them    
    untrials = 1;
    
    % instrument to give us trial and example counts
    trialCounts = zeros(nconds,1);
    imgCounts   = zeros(nconds,1);
    
    for nt=1:1:ntrials
            
      len   = info(nt).len;
      cond  = info(nt).cond;
      tinfo = info(nt);
      tinfo.len = 1; % so that it can be used for each image

      if cond > 0
	s = trialBegin(cond); e = trialEnd(cond);
	trialCounts(cond) = trialCounts(cond) + 1;
	imgCounts(cond)   = imgCounts(cond) + (e-s+1);
      else
	s = 1; e = len;
      end

      for t=s:1:e
	udata{untrials} = data{nt}(t,:);
	uinfo(untrials) = tinfo;
	untrials = untrials + 1;
      end
    end
    data = udata;
    info = uinfo;
    meta.ntrials = untrials-1;

    for c=1:nconds
      fprintf('condition %d\t%d trials\t%d images\n',c,trialCounts(c),imgCounts(c));
    end
    
    % info+data+meta were transformed into a pseudo-dataset
    % with as many trials as images (minus the ones excluded)
    % and where each trial has length 1
    
    % Now replace each trial/image by an average of itself and its
    % neighbours
    ndata = cell(size(data));
    
    for t=1:meta.ntrials
      ndata{t} = zeros(size(data{t}));

      wbegin = t-radius;
      if wbegin<=0; wbegin = t; end
      wend   = t+radius;
      if wend>meta.ntrials; wend = t; end
    
      for k=wbegin:1:wend; ndata{t} = ndata{t} + data{k}; end
      ndata{t} = ndata{t} / (wend-wbegin+1);
    end
      
    data = ndata;
          
   otherwise
    % set condition used to determine data abstraction
    condition = 2; % fix this for a while
    info = rinfo{condition}; data = rdata{condition}; meta = rmeta{condition};
  end


% Apply preprocessing steps to IDM  
%
% The wrapper file contains a cell array specifying a sequence of
% transformations (one per cell). Each transformation comprises:
% 1) function to apply to IDM
% 2) arguments required in addition to IDM
%
% This function goes through the cell array and applies each transformation

function [info,data,meta] = transformIDM_preProcess(info,data,meta,processSteps)
  nSteps = length(processSteps);

  for s=1:nSteps

    step = processSteps{s};
    stepName = step{1};
    stepArgs = step{2};
    nArgs    = length(stepArgs);

    fprintf('\tloadSubjectDataMult: transformIDM_preProcess: %s\n',stepName);
    
    % Create temporary variables to hold arguments
    % Allows us to create an argument string containing only
    % variable names
    argString = 'info,data,meta'; % these are always needed
    for a=1:nArgs
      eval(sprintf('tmpVar%d = stepArgs{a};',a));
      argString = [argString,sprintf(',tmpVar%d',a)];
    end

    % Now create a command for the preprocessing step using the argString
    % and execute it
    cmd = sprintf('[info,data,meta]=%s(%s);',stepName,argString);
    eval(cmd);

    if 0
      % uncomment to have printouts of what is being called in each
      % step, and with what arguments      
      for a=1:nArgs
	if isnumeric(stepArgs{a})
	  isfloat = stepArgs{a} - round(stepArgs{a});
	  if isfloat
	    fprintf('\t%f\n',stepArgs{a});
	  else
	    fprintf('\t%d\n',stepArgs{a});
	  end	
	elseif ischar(stepArgs{a})
	  fprintf('\t%s\n',stepArgs{a});
	else
	  fprintf('\tCannot handle type of argument %d\n',a);
	end
      end
    end
  
  end

