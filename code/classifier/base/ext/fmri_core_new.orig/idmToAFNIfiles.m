function idmToAFNIfiles(varargin)

%Purpose:
%
%   Convert the dataset(in the format of info, data and meta) into AFNI's format.
%   This function is at first designed to store the functional data into AFNI's format.
%   If the user does specify the snapshot number, the chosen snapshot will be 
%   stored into a functional dataset, otherwise they will be stored into an anatomical
%   dataset, in order to display the time series curve using AFNI. If you have a 
%   functional matrix already, you can use WriteMatrix directly to write it into AFNI's 
%   format(for details, see README.convert).
%
%     
%Input Parameters:
%
%   anatDataset: the name of the existed anatomical dataset of this subject
%   coregHeader: the coregistrated head file for this subject
%   outputDataset: the prefix name of the new functional dataset
%   trialNum: trial number
%   snapshotNum: snapshot number
%
%
%Output Parameters:
%
%   None
%
%
%Examples:
%
%   [i,d,m] = loadSubjectdata('data-categories', '05886', {'LIT' 'CALC'},1);
%   idmToAFNIfiles(i, d, m, '05886_oblique_spgr+orig','05886+orig.HEAD','New', 3);
%   create the AFNI files to display time series of trial 3. Store the
%   results on the current directory as files beginning with prefix 'New'.
%
%More Info:
%   idmToAFNIfiles(i,d,m,anatDataset,coregHeader,outputDataset,trialNum,<snapshotNum>)
%   Create the files needed to display the data using the AFNI program.
%   You should add the AFNI directory into your path to let this function work well.
%   This function will call the functions in the subdirectory named 'Convert'.
%   For details, see /afs/cs/project/theo-72/fmri_core/Convert/README.convert 
%
%   anatDataset: the AFNI anatomical dataset for the subject of interest. 
%   This must be on the current directory. If you haven't this dataset,
%   see Hint 3 in README.convert
%   coregHeader: name of the coregistrated head file for this subject.
%   If you haven't this headerfile, first try to get it in
%   /afs/cs/project/theo-70/Src/headers/. If failed, see Hint 4 in README.convert.
%   outputDataset: name of the output dataset (prefix for two new
%   files that will be created on the same directory)
%   trialNum: number of the trial to be displayed (if you want to display
%   multiple trials, first use transformIDM_selectTrials and
%   transformIDM_mergeTrials)
%   snapshotNum: number of the snapshot within the trial to be displayed.
%   This is optional.  If not given, the entire trial will be displayed.
%  
%	
%See Also:
%
%   BrikInfo, WriteBrik, BrikRead
%
%
%Author:
%
%   Xuerui Wang(xuerui@cs.cmu.edu)
%
%
%Date:
%
%   Created on Sat Sep. 4, 2002
%   Updated on Sun Sep. 7, 2002
%   Updated on Tue Sep. 17, 2002
% - 18 Sep 2002 - fp - vectorised the code that creates the 3D matrix.
%                      Still have to do the 4D one.
% - 19 Sep 2002 - fp - added image text tag support
  
snapshotNum = 0;

l = length(varargin);
if l < 8
    fprintf(1,'syntax: idmToAFNIfiles(i, d, m, anatDataset, outputDataset, trialNum, <snapshotNum>)\n\n');
    fprintf(1,'- anatDataset: the name of the existed anatomical dataset\n');
    fprintf(1,'coregHeader: the coregistrated head file for this subject\n')
    fprintf(1,'- outputDataset: the prefix name of the new functional dataset\n');
    fprintf(1,'- trialNum: trial number\n');
    fprintf(1,'- snapshotNum: snapshot number\n');
    fprintf(1,'\n');
    return;
else
    info = varargin{1};
    data = varargin{2};
    meta = varargin{3};
    AnatDatasetName = varargin{4};
    coregHeader     = varargin{5};
    NewPrefix       = varargin{6};
    trialNum        = varargin{7};
    snapshotNum     = varargin{8};
    doTagging = 0;
    if l > 8
      labels    = varargin{9};
      doTagging = 1;
    end
end

timepoints = size(data{trialNum},1);
nvoxels = meta.nvoxels;

if(snapshotNum ~= 0)
  % code as of sep 18 - a bit slow because of the loop over all
  % voxels and does not work with all studies
  % (remember that meta has the dimensions of the grid the data came from)
  %
  %    M = zeros(64,64,16);
  %
  %    for(i = 1:1:nvoxels)
  %        j = meta.colToCoord(i,1);
  %        k = meta.colToCoord(i,2);
  %        m = meta.colToCoord(i,3);
  %        M(j,k,m) = data{trialNum}(i);
  %    end
  
  % vectorised version
  dimx = meta.dimx; dimy = meta.dimy; dimz = meta.dimz; timepoints = 1;
  M = zeros(dimx,dimy,dimz);
  % from coordinates, compute indices in a vector version
  % of a 3D matrix (i.e all columns of the matrix lined into a vector);
  xcoords = meta.colToCoord(:,1);
  ycoords = meta.colToCoord(:,2);
  zcoords = meta.colToCoord(:,3);
  indices = sub2ind([dimx dimy dimz],xcoords,ycoords,zcoords);
  % and store in the matrix
  M(indices) = data{trialNum}(snapshotNum,:);
  
else
  dimx = meta.dimx; dimy = meta.dimy; dimz = meta.dimz; timepoints = 1;
  M = zeros(dimx,dimy,dimz,timepoints);
  
  for(i = 1:1:nvoxels)
    j = meta.colToCoord(i,1);
    k = meta.colToCoord(i,2);
    m = meta.colToCoord(i,3);
    for(n =1: 1: timepoints)
      M(j,k,m,n) = data{trialNum}(n,i);
    end
  end
end

% Tagging

if doTagging

  ntags = length(labels);
  zpossible = flipud(sort(unique(zcoords))); % highest z first
  colorNumber = 0.5; 
  
  % for now just tag the first time point  
  for t=1:1
    
    for z = 1:min(ntags,dimz)
      tags = labels{z}; l = length(tags);
      if l > 1; doBelow=1; else doBelow=0; end
      
      pictograms = makeTag(tags{1});
      
      pl = length(pictograms);
      xidx = 1; yidx = 1; minzidx = min(zcoords); maxzidx = max(zcoords);
      for p=1:pl
	M(xidx:xidx+2,yidx:yidx+4,zpossible(z),t) = pictograms{p}'*colorNumber;
	xidx = xidx + 4;
	if xidx > dimx; break; end
      end
    
      if doBelow
	pictograms = makeTag(tags{2});
	
	pl = length(pictograms);
	xidx = 1 ; yidx = meta.dimy-5; minzidx = min(zcoords); maxzidx = max(zcoords);
	for p=1:pl
	  M(xidx:xidx+2,yidx:yidx+4,zpossible(z),t) = pictograms{p}'*colorNumber;
	  xidx = xidx + 4;
	  if xidx > dimx; break; end
	end
      end
  
    end; % for over slices
  end; % for over time points
    
end

%
% Create the AFNI files
%

WriteMatrix(M, AnatDatasetName, NewPrefix);
cmd = sprintf('!cp %s %s+orig.HEAD',coregHeader,NewPrefix);
eval(cmd);
NewDataset = sprintf('%s+orig',NewPrefix);
ModifyHead(NewDataset,M);
cmd = sprintf('!3dnewid %s',NewDataset);
eval(cmd);

return;



% only upper case or decimal
  
function [pictograms] = makeTag(tagText)
  
  tagText = upper(tagText);
  l = length(tagText);
  for p=1:l
    c = tagText(p);  
    switch c
     case {'A'}; pc = [[1 1 1];[1 0 1];[1 1 1];[1 0 1];[1 0 1]];
     case {'B'}; pc = [[1 1 0];[1 0 1];[1 1 0];[1 0 1];[1 1 0]];
     case {'C'}; pc = [[1 1 1];[1 0 0];[1 0 0];[1 0 0];[1 1 1]];
     case {'D'}; pc = [[1 1 0];[1 0 1];[1 0 1];[1 0 1];[1 1 0]];
     case {'E'}; pc = [[1 1 1];[1 0 0];[1 1 1];[1 0 0];[1 1 1]];
     case {'F'}; pc = [[1 1 1];[1 0 0];[1 1 0];[1 0 0];[1 0 0]];
     case {'G'}; pc = [[1 1 1];[1 0 0];[1 1 1];[1 0 1];[1 1 1]];     
     case {'H'}; pc = [[1 0 1];[1 0 1];[1 1 1];[1 0 1];[1 0 1]];
     case {'I'}; pc = [[0 1 0];[0 1 0];[0 1 0];[0 1 0];[0 1 0]];
     case {'J'}; pc = [[0 1 0];[0 1 0];[0 1 0];[0 1 0];[1 1 0]];
     case {'K'}; pc = [[1 0 1];[1 1 0];[1 0 0];[1 1 0];[1 0 1]];
     case {'L'}; pc = [[0 1 0];[0 1 0];[0 1 0];[0 1 0];[0 1 1]];
     case {'M'}; pc = [[1 0 1];[1 1 1];[1 0 1];[1 0 1];[1 0 1]];
     case {'N'}; pc = [[1 0 1];[1 0 1];[1 1 1];[1 0 1];[1 0 1]];
     case {'O'}; pc = [[1 1 1];[1 0 1];[1 0 1];[1 0 1];[1 1 1]];
     case {'P'}; pc = [[1 1 0];[1 0 1];[1 1 0];[1 0 0];[1 0 0]];
     case {'Q'}; pc = [[1 1 1];[1 0 1];[1 0 1];[1 0 1];[1 1 1]];
     case {'R'}; pc = [[1 1 0];[1 0 1];[1 1 0];[1 1 0];[1 0 1]];
     case {'S'}; pc = [[1 1 1];[1 0 0];[1 1 1];[0 0 1];[1 1 1]];
     case {'T'}; pc = [[1 1 1];[0 1 0];[0 1 0];[0 1 0];[0 1 0]];
     case {'U'}; pc = [[1 0 1];[1 0 1];[1 0 1];[1 0 1];[1 1 1]];
     case {'V'}; pc = [[1 0 1];[1 0 1];[1 0 1];[1 0 1];[0 1 0]];
     case {'X'}; pc = [[1 0 1];[1 0 1];[0 1 0];[1 0 1];[1 0 1]];
     case {'Y'}; pc = [[1 0 1];[1 0 1];[0 1 0];[0 1 0];[0 1 0]];
     case {'W'}; pc = [[1 0 1];[1 0 1];[1 0 1];[1 0 1];[0 1 0]];
     case {'Z'}; pc = [[1 1 1];[0 0 1];[0 1 0];[1 0 0];[1 1 1]];
     case {'0'}; pc = [[1 1 1];[1 0 1];[1 0 1];[1 0 1];[1 1 1]];
     case {'1'}; pc = [[0 1 0];[0 1 0];[0 1 0];[0 1 0];[0 1 0]];
     case {'2'}; pc = [[1 1 1];[0 0 1];[1 1 1];[1 0 0];[1 1 1]];
     case {'3'}; pc = [[1 1 1];[0 0 1];[1 1 1];[0 0 1];[1 1 1]];
     case {'4'}; pc = [[1 0 1];[1 0 1];[1 1 1];[0 0 1];[0 0 1]];
     case {'5'}; pc = [[1 1 1];[1 0 0];[1 1 1];[0 0 1];[1 1 1]];
     case {'6'}; pc = [[1 0 0];[1 0 0];[1 1 1];[1 0 1];[1 1 1]];
     case {'7'}; pc = [[1 1 0];[0 1 0];[0 1 0];[0 1 0];[0 1 0]];
     case {'8'}; pc = [[1 1 1];[1 0 1];[1 1 1];[1 0 1];[1 1 1]];
     case {'9'}; pc = [[1 1 1];[1 0 1];[1 1 1];[0 0 1];[0 0 1]];      
     case {'_',' '}; pc = [[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0]];
     case {'-'}; pc = [[0 0 0];[0 0 0];[1 1 1];[0 0 0];[0 0 0]];
     otherwise
      fprintf('error: makeTag cannot deal with character %c, replacing by space\n',c);
      pc = [[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0]];
    end
    
    pictograms{p} = [pc];
  end
    