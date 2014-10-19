% TESTED 8/20/02.
%
% working example illustrating animation and other data visualization
%
% Note to execute this you must be on a machine containing the data
% for the 'data-categories' study for subject 05886. 


study   = 'data-categories';
subject = '05886'  
rois = {'LIT' 'CALC'};

% if you want see really beautiful pictures, then load (more slowly) all
% of the ROIs for this subject:
% rois = {'ACING' 'LDLPFC' 'LIPL' 'LOPER' 'LSGA' 'LTRIA' 'RDLPFC' 'RIPL' ...
%	'ROPER' 'RSGA' 'RTRIA' 'CALC' 'LFEF' 'LIPS' 'LPPREC' 'LSPL' 'OP' ...
%	'RFEF' 'RIPS' 'RPPREC' 'RSPL' 'SMA' 'LIES' 'LIT' 'LSES' 'LT' 'RIES' ...
%	'RIT' 'RSES' 'RT' 'SMFP'}; 


% load the data.  returns cell arrays of IDMs
[i,d,m] = loadSubjectdata(study,subject,rois);

% plot an image showing the brain activity for the 7th snapshot in trial 3,
% displaying the z=10 slice
n=7; trialNum=3; z=10; 
plotSnapshot(i,d,m,trialNum,z,n,0,0);

% plot the entire time series for trial 3, for each voxel in each brain slice
plotTrial(i,d,m,3)

% watch a movie of the brain activity for brain slice z=10 of trial 3
M=animateTrial(i,d,m,3,10);

% watch it again
movie(M);

% watch a movie of the brain activity for six z slices, begining with z=5, of
% trial 3
M=animate6Trial(i,d,m,3,5);


% watch a movie of the brain activity for all 16 z slices of
% trial 3
M=animate16Trial(i,d,m,3);

% erase the big movie data structure
clear M;
