% Merge two datasets (info+data+meta)
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
% - 24 Oct, 2002 Xuerui - kept study, subject and other fields of meta, too.


function [info,data,meta] =mri_mergeIDM( info1,data1,meta1,info2,data2,meta2 )

% test that # trials, lengths and conditions match
  if meta1.ntrials ~= meta2.ntrials
    fprintf(1,'mri_compareDatasets: error: # trials does not match.\n');
    result = 0; return;
  end
  
  meta = meta1;
  ntrials      = meta1.ntrials;
  meta.ntrials = ntrials;
  
  for nt=1:1:ntrials  
    if info1(nt).len ~= info2(nt).len
      fprintf(1,'mri_compareDatasets: error: info len mismatch in trial %d\n',nt);
      result = 0; return;
    end
    if info1(nt).cond ~= info2(nt).cond
      fprintf(1,'mri_compareDatasets: error: info cond mismatch in trial %d\n',nt);
      result = 0; return;
    end      
  end

  % a few parameters and allocations
  
  info = info1; % info1=info2 must match
  data = cell(ntrials,1);
  nvoxels1 = size(data1{1},2);
  nvoxels2 = size(data2{1},2);
  nvoxels  = nvoxels1 + nvoxels2;
  
  meta.nvoxels = nvoxels;
  meta.dimx = meta1.dimx;
  meta.dimy = meta1.dimy;
  meta.dimz = meta1.dimz;
  
  % a dirty trick, but helps in merging 1 voxel average ROIs
  if nvoxels1 == 1 & nvoxels2 == 1
    for nt=1:1:ntrials
%      data{nt} = zeros(info(nt).len,2);
      
    end
  end
  
  % merge data
  % For each trial
  % - create a temporary array capable of containing the transposes
  %   of the data for both trials, in addition to three columns with
  %   their coordinates
  % - sort them by xyz
  % - store the remainder in the result
    
  for nt=1:1:ntrials
    data{nt} = zeros(nvoxels,info(nt).len+3);
    tlen = info(nt).len;
    % copy transpose of trial data into destination
    % i.e. T1 T2 end up as
    % T1'
    % T2'   
    data{nt}(1:1:nvoxels1,1:1:tlen) = data1{nt}';
    data{nt}(nvoxels1+1:1:nvoxels1+nvoxels2,1:1:tlen) = data2{nt}'; 
    % copy coordinate information to last three columns
    % T1'C1
    % T2'C2
    data{nt}(1:1:nvoxels1,tlen+1:1:tlen+3) = meta1.colToCoord(1:1:nvoxels1,:);
    data{nt}(nvoxels1+1:1:nvoxels1+nvoxels2,tlen+1:1:tlen+3) = meta2.colToCoord(1:1:nvoxels2,:);
    % sort by zyx
    data{nt} = sortrows(data{nt},[tlen+3 tlen+2 tlen+1]);    
    % data{nt} = data{nt}';
    data{nt} = (data{nt}(:,1:1:tlen))';
  end
  
  % create appropriate mappings
  meta.colToCoord = zeros(nvoxels1+nvoxels2,3);
  meta.colToCoord(1:1:nvoxels1,:) = meta1.colToCoord(1:1:nvoxels1,:);
  meta.colToCoord(nvoxels1+1:1:nvoxels1+nvoxels2,:) = meta2.colToCoord(1:1:nvoxels2,:);
  meta.colToCoord = sortrows(meta.colToCoord,[3 2 1]);
  meta.coordToCol = zeros(meta.dimx,meta.dimy,meta.dimz);

  for col=1:1:meta.nvoxels
    cx = meta.colToCoord(col,1);
    cy = meta.colToCoord(col,2);
    cz = meta.colToCoord(col,3);    
    meta.coordToCol(cx,cy,cz) = col;
  end
    

  
  
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
