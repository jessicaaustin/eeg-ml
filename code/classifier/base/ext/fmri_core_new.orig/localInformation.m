% Fetches information about:
% - location of code/data (or links to them) in this machine
% - which studies/subjects are available
%
% The script can be configured for someone's particular machine,
% which is useful if the person happens to have a partial copy of
% the dataset. Check the code for more information on how to
% do this, as there are sections for different people in there already.
%
% In order to do this, pick a place in your machine for storing all
% fmri data and code (let's call it <LOCAL PATH>). Then create
% two subdirectories in it:
% 1) fmri_core
% 2) fmri_data
%
% and copy
% - the /afs/cs/project/theo-72/fmri_core contents to 1)
% - any dataset directories (e.g. data-starplus-sp,
%   data-categories) to 2)
%
% Then set the machine specific section below appropriately.
%
%
%
% In:
% - dataset - a dataset name, of the form "data-<study name>"
%
% Out:
% - information - a structure containing the information,
%
%  What host are we in
%  - information.host = HOSTNAME;

%  What directory is data in on this HOSTNAME
%  - information.dataplace = dataplace;

%  What directory is code in on this HOSTNAME
%  - information.codeplace = codeplace;

%  What subjects are available on this machine?
%  - information.subjects  = subjects;

%  Which ROIs are available?
%  - information.rois      = rois;

%  What are the 3D dimensions of this dataset (returns [x y z])?
%  - information.dimensions = dimensions;

%  How many conditions, including fixation?
%  - information.nconds    = nconds;
% 
% Example
%
% information = localInformation('data-starplus-sp');
%
% History:
%
% - 17 Mar 2005 - fpereira - adapted from old mri_reference.m

function [information] = localInformation(varargin)
  
  l = length(varargin);
  if l < 1; fprintf('syntax: localInformation(<study name>)\n'); return;
  else
    dataset = varargin{1};
  end
 
  
  %% What machine are we running on ?

  global HOSTNAME;
  if isempty(HOSTNAME)
    % can set this globally (recommended if running with jvm)
    % or ask operating system (works fine when 'matlab -nojvm')
    [s,HOSTNAME] = unix('hostname');
    HOSTNAME = lower(deblank(HOSTNAME));
  end;

  %%
  %% Machine specific settings
  %%
  
  %% locations for code/data that depend on which 
  %% machine you are running this function on on)
  
  % numeric flag used to denote subjects/studies
  reduce    = 0; % all subjects available for all studies
  codeplace = '/afs/cs/project/theo-72/fmri_core_new';
  
  switch HOSTNAME
   case {'human','human.learning.cs.cmu.edu'}
    % Tom's laptop
    dataplace = 'E:\fMRI\Data';
    codeplace = '';
    reduce = 3; % laptop specific reduced set of subjects/ROIs    
   case {'animal','animal.learning.cs.cmu.edu'}
    % Tom's animal laptop
    dataplace = 'E:/fMRI/data';
    codeplace = 'C:\Documents and Settings\mitchell\My Documents\RESEARCH\fMRI\fmri_core_new';
    reduce = 3; % laptop specific reduced set of subjects/ROIs    
   case {'natural','natural.learning.cs.cmu.edu'}
    % Tom's laptop in Windows2000
    dataplace = 'F:\fMRI_data';
    codeplace = '';
    reduce = 3; % laptop specific reduced set of subjects/ROIs    
   case {'pingo','pingo.pc.cs.cmu.edu'}
    % francisco's laptop
    dataplace = '/home/fpereira/projecto/FMRI/fmri_data';
    codeplace = '/home/fpereira/projecto/FMRI/fmri_core_new';
    reduce = 2; % laptop specific reduced set of subjects/ROIs
   case {'gs2011','gs2011.sp.cs.cmu.edu'}
    % rebecca's machine
    dataplace = '/usr0/rah/FMRIdata';
    codeplace = '';
    reduce = 1; % workstation specific reduced set of subjects/ROIs
    % add a comment
   case {'gs233','gs233.sp.cs.cmu.edu'}
    % stefan's machine
    dataplace = '/usr2/proiecte/fmri/fmri_data';
    codeplace = '';
    reduce = 3; % workstation specific reduced set of subjects/ROIs
    % add a comment
    % hope it works :)
   case {'gs3098', 'gs3098.sp.cs.cmu.edu'}
    % indra's machine
    dataplace = '/usr0/fmri/fmri_data';
    codeplace = '';
   case {'arjuna', 'arjuna.wv.cc.cmu.edu'}
    % indra's laptop
    dataplace = '/Users/irustand/fmri_data';
    codeplace = '';
   case {'wew.fmri.cs.cmu.edu'}
    % Wei's desktop
    dataplace = '/usr1/fMRI/Data';
    codeplace = '/usr1/fMRI/fmri_core_new';
    reduce = 3; % laptop specific reduced set of subjects/ROIs
   
   case{'node2'}
    dataplace='/scratch';
    codeplace='/usr/cluster/projects/cat_scripts/fmri_core_new';
   case{'node6'}
    dataplace='/scratch';
    codeplace='/usr/cluster/projects/cat_scripts/fmri_core_new';
   case{'node5'}
    dataplace='/scratch';
    codeplace='/usr/cluster/projects/cat_scripts/fmri_core_new';
     
       case{'node7'}
    dataplace='/scratch';
    codeplace='/usr/cluster/projects/cat_scripts/fmri_core_new';
       case{'node3'}
    dataplace='/scratch';
    codeplace='/usr/cluster/projects/cat_scripts/fmri_core_new';
   case{'node10'}
    dataplace='/scratch';
    codeplace='/usr/cluster/projects/cat_scripts/fmri_core_new';
     case{'node11'}
    dataplace='/scratch';
    codeplace='/usr/cluster/projects/cat_scripts/fmri_core_new';
      case{'node14'}
    dataplace='/scratch';
    codeplace='/usr/cluster/projects/cat_scripts/fmri_core_new';
   case{'node15'}
    dataplace='/scratch';
    codeplace='/usr/cluster/projects/cat_scripts/fmri_core_new';
   case{'node4'}
    dataplace='/scratch';
    codeplace='/usr/cluster/projects/cat_scripts/fmri_core_new';
   case{'node13'}
    dataplace='/scratch';
    codeplace='/usr/cluster/projects/cat_scripts/fmri_core_new';
     case{'node9'}
    dataplace='/scratch';
    codeplace='/usr/cluster/projects/cat_scripts/fmri_core_new';
   case{'node8'}
    dataplace='/scratch';
    codeplace='/usr/cluster/projects/cat_scripts/fmri_core_new';
   case{'node12'}
    dataplace='/scratch';
    codeplace='/usr/cluster/projects/cat_scripts/fmri_core_new';
   case{'node16'}
    dataplace='/scratch';
    codeplace='/usr/cluster/projects/cat_scripts/fmri_core_new';
   case{'node16'}
    dataplace='/scratch';
    codeplace='/usr/cluster/projects/cat_scripts/fmri_core_new';
   case {'multivac.cald.cs.cmu.edu'}

    % MULTIVAC - main data store
    
    switch dataset
     case {'data-brainlex','data-brainlexNew','data-brainlexNotNormalized','data-brainlex-org'...
      'data-categories','data-exemplar','data-goldcategories','data-sixcategories'...
      'data-starplus-ps','data-starplus-sp','data-starplus-spps','data-syntamb2'...
      'data-threecategories','data-threecategories-old','data-twocategories','data-fourcategories',...
      'data-biEng', 'data-biPort'}
      dataplace = '/usr3/fmri/fmri_data';
     case {'justTheFactsMam'}
      dataplace = '/usr3/fmri/fmri_data';
    end
    
   otherwise

    fprintf('localInformation: ERROR:\n');
    fprintf(' machine %s is unknown, please add a section for it in localInformation.m\n',HOSTNAME);
    pause;return;
    
  end

  
  %
  % Now specify which subjects are available in each dataset in your machine,
  % by setting reduce to a number > 0. Leave unchanged if you have
  % all.
  % 
  % NOTE: undefined if NO subjects exist on the machine. Not really
  % a problem because the person shouldn't be calling this if they
  % don't have the dataset on their machine. But it should be fixed
  % somehow...
  %
  % information for each dataset
  % - subjects   - a cell array with subject IDs
  % - ROIs       - a cell array with ROI IDs
  % - dimensions - the dimensions of the 3D volume the data is in
  % - nconds     - the # of conditions, including fixation
  %
    
  % subjects available
  switch dataset

   case 'data-biEng'  
    subjects = {'01701B' '01708B' '01723B' '01730B' '01751B' '01765B' '01771B' '01776B' '01778B' '01783B'};
    switch reduce
     case 1
      subjects = {'01701B'};
     otherwise
      subjects = sort(subjects);
    end       

    % ROIS not available
    rois = {'Y' 'VPSC' 'MPSC'};
    dimensions = [64 64 16]; nconds = 2;
      
   case 'data-biPort'  
    subjects = {'01701B' '01708B' '01723B' '01730B' '01751B' '01765B' '01771B' '01776B' '01778B' '01783B'};
    switch reduce
     case 1
      subjects = {'01701B'};
     otherwise
      subjects = sort(subjects);
    end       

    % ROIS not available
    rois = {'Y' 'VPSC' 'MPSC'};
    dimensions = [64 64 16]; nconds = 2;
      
   case 'data-twocategories'  
    subjects = {'01471B' '01480B' '01481B' '01482B' '01485B' '01492B' '01544B'};
    switch reduce
     case 1
      subjects = {'01471B'};
     otherwise
      subjects = sort(subjects);
    end       

    % ROIS not available
    rois = {'Y' 'VPSC' 'MPSC' 'BRAIN_FDS'};
    dimensions = [64 64 16]; nconds = 2;
      
   case 'data-fourcategories'  
    subjects = {'01244B' '01269B' '01288B' '01286B' '01281B' '01306B' ...
        '01310B' '01311B' '01318B' '01326B' '01362B' '01357B' '01370B' '01375B'};
    switch reduce
     case 1
      subjects = {'01244B'};
     otherwise
      subjects = sort(subjects);
    end       

    % ROIS not available
    rois = {'BRAIN_FDS'};
    dimensions = [64 64 16]; nconds = 2;
    
   case 'data-starplus-sp'  
    subjects_with_baddata={ '04772_20_sp3terC'};
    subjects = {'04847_20_sp3terC' '04805_30_sp3terC' ...
		'05005_30_sp3terC' '05093_20_sp3terC' '04799_20_sp3terC' ...
		'04820_20_sp3terC' '04958_30_sp3terC' '05018_30_sp3terC' ...
		'05393_30_sp3terC' '05099_20_sp3terC' '05643_20_sp3terC' ...
		'05675_10_sp3terC' '05680_10_sp3terC' '05695_10_sp3terC' ...
		'05710_20_sp3terC' '05131_30_sp3terC'};
 
    switch reduce
     case 1
      subjects = {'04847_20_sp3terC'};
     case 2
      subjects = {'04847_20_sp3terC','04799_20_sp3terC'};
     case 3
      subjects = {'04847_20_sp3terC','04820_20_sp3terC','04799_20_sp3terC'};
     otherwise
      subjects = sort(subjects);
    end

    % ROIS available
    rois = {'CALC' 'LFEF' 'LIPL' 'LIT' 'LPPREC' 'LSPL' 'LTRIA' 'RFEF' 'RIPS' 'ROPER' 'RSGA' 'RT' 'SMA' 'LDLPFC' 'LIPS' 'LOPER' 'LSGA' 'LT' 'RDLPFC' 'RIPL' 'RIT' 'RPPREC' 'RSPL' 'RTRIA'};

    dimensions = [64 64 8]; nconds = 3;
   
   case 'data-starplus-ps'
    subjects = {'04772_30_ps3terC' '04799_30_ps3terC' '04805_20_ps3terC' ...
		'04820_30_ps3terC' '04847_30_ps3terC' '04958_20_ps3terC' ...
		'05005_20_ps3terC' '05018_20_ps3terC' '05093_30_ps3terC' ...
		'05099_30_ps3terC' '05131_20_ps3terC' '05393_20_ps3terC' ...
		'05643_10_ps3terC' '05675_20_ps3terC' '05680_20_ps3terC' ...
		'05695_20_ps3terC' '05710_10_ps3terC'};

    switch reduce
     case 1
      subjects = {'04847_30_ps3terC'};
     case 2
      subjects = {'04847_30_ps3terC','04799_30_ps3terC'};
     case 3
      subjects = {'04847_30_ps3terC','04820_30_ps3terC','04799_30_ps3terC'};
     otherwise
      subjects = sort(subjects);
    end

    % ROIS available
    rois = {'CALC' 'LFEF' 'LIPL' 'LIT' 'LPPREC' 'LSPL' 'LTRIA' 'RFEF' 'RIPS' 'ROPER' 'RSGA' 'RT' 'SMA' 'LDLPFC' 'LIPS' 'LOPER' 'LSGA' 'LT' 'RDLPFC' 'RIPL' 'RIT' 'RPPREC' 'RSPL' 'RTRIA'};
    
   dimensions = [64 64 8]; nconds = 3;

    case 'data-starplus-spps'
    subjects = {'04772' '04799' '04805' ...
		'04820' '04847' '04958' ...
		'05005' '05018' '05093' ...
		'05099' '05131' '05393' ...
		'05643' '05675' '05680' ...
		'05695' '05710'};

    switch reduce
     case 1
      subjects = {'04847'};
     case 2
      subjects = {'04847','04799'};
     case 3
      subjects = {'04847','04820','04799'};
     otherwise
      subjects = sort(subjects);
    end

    % ROIS available
    rois = {'CALC' 'LFEF' 'LIPL' 'LIT' 'LPPREC' 'LSPL' 'LTRIA' 'RFEF' 'RIPS' 'ROPER' 'RSGA' 'RT' 'SMA' 'LDLPFC' 'LIPS' 'LOPER' 'LSGA' 'LT' 'RDLPFC' 'RIPL' 'RIT' 'RPPREC' 'RSPL' 'RTRIA'};
    
   dimensions = [64 64 8]; nconds = 3; 
   % loadSubjectdata, please check loadSubjectdataSPPS(Subjects,ROIs) in README.Function
   
   case 'data-syntamb2'
    subjects = {'02945_22' '02953_32' '02964_22' '02999_24' '03002_31' '03003_32'};

    switch reduce
     case 1
      subjects = {'02945_22'};
     case 4
      subjects = {'02945_22' '03003_32'};
     otherwise
      subjects = sort(subjects);
    end

    rois = {'LB' 'LT' 'RB' 'RT'};
    dimensions = [64 64 7]; nconds = 7;
    
    
   case {'data-brainlex','data-brainlexNew','data-brainlexNotNormalized'}

    % Data was extracted from
    % binaries that were corrected for timing of acquisition of
    % slices and had missing images smoothed over (6 subjects total).
    
    % In addition, notNormalized is not normalized % above fixation.
    % For some reason, the detrend.dat binary for 08330 was not
    % correcte or smoothed, so that was not extracted.
    
    subjects = {'08057' '08170' '08179' '08240' '08276' '08298'};    
    
    switch reduce
     case 1
      subjects = {'08057'};
     case 2
      subjects = {'08057' '08298'};
     otherwise
      subjects = sort(subjects);
    end
    
    % minimal set of ROIs that all the subjects have
    rois = {'ACING' 'LIES' 'LIT' 'LSES' 'LT' 'RIES' 'RIT' 'RSES' ...
	    'RT' 'SMFP' 'LDLPFC' 'LIPL' 'LOPER' 'LSGA' 'LTRIA' ...
	    'RDLPFC' 'RIPL' 'ROPER' 'RSGA' 'CALC' 'LFEF' 'LIPS' ...
	    'LPPREC' 'LSPL' 'OP' 'RFEF' 'RIPS' 'RSPL' 'SMA' 'RTRIA' ...
	    'LCBELL' 'RCBELL' 'RPPREC'};

    % There's also ALLBRAIN, with more voxels than all of these put together.
    % rois = {rois{:},'ALLBRAIN'};
    % and the cerebellum, if you really want it
    % rois = {rois{:},'LCBELL','RCBELL'};

    dimensions = [64 64 14]; nconds = 4;

    %
    % Block design category studies 
    % (12, 6 and 3 categories with 1,  2 and 4 blocks of stimuli of each category,respectively)
    %

   case {'data-categories','data-sixcategories','data-threecategories'}

    rois = {'ACING' 'CALC' 'OP' 'SMA' 'SMFP' 'LIES' 'RIES' 'LIT' 'RIT' ...
	    'LT' 'RT' 'LIPL' 'RIPL' 'LOPER' 'ROPER' 'LTRIA' 'RTRIA' ...
	    'LSGA' 'RSGA' 'LFEF' 'RFEF' 'LPPREC' 'RPPREC' 'LDLPFC' 'RDLPFC' ...
	    'LIPS' 'RIPS' 'LSES' 'RSES' 'LSPL' 'RSPL'};

    % There's also ALLBRAIN, with more voxels than all of these put together.
    % rois = {rois{:},'ALLBRAIN'};
    % and the cerebellum, if you really want it
    % rois = {rois{:},'LCBELL','RCBELL'};
    
    switch dataset
      
     case {'data-categories'}

      subjectsWithBadData = {'05688'};
      subjects = {'05483' '05502' '05549' '05675' '05695' '05797' '05813' '05833' '05838' '05868' '05874' '05886'};

      switch reduce
       case 1
	subjects = {'05886'};
       case 2
	subjects = {'05868' '05688' '05833'};
       otherwise
	subjects = sort(subjects);
      end
      dimensions = [64 64 16]; nconds = 13;
      
     case {'data-sixcategories'}

      subjects = {'233B' '86B' '77B' '329B' '332B' '424B' '474B' '496B'};
    
      switch reduce
       otherwise
	subjects = sort(subjects);
      end
      dimensions = [64 64 16]; nconds = 7;
       
     case {'data-threecategories'}

      subjects = {'354B','357B','362B','367B','371B'};

      switch reduce
       otherwise
	subjects = sort(subjects);
      end
      dimensions = [64 64 16]; nconds = 4;

    end

   case {'justTheFactsMam'}
    
    % used for calls that don't want study information, just
    % information about paths in this hostname;
    subjects = {}; rois = {};
    dimensions = []; nconds = 0;
    
   otherwise
    
  end

  information.host = HOSTNAME;
  information.dataplace  = dataplace;
  information.codeplace  = codeplace;
  information.subjects   = subjects;
  information.rois       = rois;
  information.dimensions = dimensions;
  information.nconds     = nconds;
