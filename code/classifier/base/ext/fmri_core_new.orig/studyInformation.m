% Fetches information about a study, namely subjects and ROIs
% available and their physical location relative to the machine
% this function is being called from.
%
% The script can be configured for someone's particular machine,
% which is useful if the person happens to have a partial copy of
% the dataset. Check the code for more information on how to
% do this.
%
% In:
% - dataset - a dataset name, of the form "data-<study name>"
%
% Out:
% - information - a structure containing the information,
%
%  What host are we in
%  - information.host = HOSTNAME;
%  What directory is this study stored on, relative to HOSTNAME
%  - information.dataplace = dataplace;
%  What subjects are available on this machine?
%  - information.subjects  = subjects;
%  Which of them should be bypassed, if any (usually, none)?
%  - information.mask      = mask;
%  Which ROIs are available?
%  - information.rois      = rois;
%  What are the 3D dimensions of this dataset (returns [x y z])?
%  - information.dimensions = dimensions;
%  How many conditions, including fixation?
%  - information.nconds    = nconds;
% 
% Dependencies:
% - none
%
% History:
%  08 Aug 02 fpereira - created core version from previous instance
%
%  14 Aug 02 tom - added global HOSTNAME, which can be set externally
%  to override the call to unix('hostname') (this unix call fails
%  under some linux implementations when using jvm)
% 
% 26 Sept 02 tom - added humanWindows as possible hostname
% 03 Oct 02 Xuerui - added the host liatris.cald.cs.cmu.edu
% 06 Oct 02 - fp - added the possibility of getting per subject ROI
% listings
% 16 Oct 02 - rah - added gs2011
% 06 Nov 02 - stefan - added gs233
% 17 Feb 03 - fp - modified lists of brainlex ROIs
% 01 Oct 03 - indra - added gs3098
% 08 Oct 03 - indra - added arjuna
% 24 Apr 04 - fp - added sections for dataset
% data-brainlexNotNormalized (check comments there for details)
% 25 Apr 04 - fp - added sections for dataset
% data-brainlexNew (check comments there for details)
% 07 Jun 04 - fp - added sections for dataset
% data-threecategoriesNew (temporary dataset while I reextract data)
% 21 Jun 04 - fp - temporary modification to extract different ROIs for data-sixcategories
% 22 Sep 04 - fp - modified the ROI list for starplus-sp and
% starplus-ps to include all the ROIs available, not just the 7
% people were using (that was just an initial suggestion). If you
% want to use only 7, please pass that list of ROIs to loadSubjectData.

%function [h,dataplace,subjects,mask,rois,dimensions,nconds] = mri_reference(dataset)
function [information] = studyInformation(varargin)
  
  l = length(varargin);
  if l < 1
    fprintf(1,'syntax: mri_reference(<dataset>,[<subject>])\n'); return;
  else
    dataset = varargin{1};
    roisPerSubject = '';
    if l > 1
      roisPerSubject = varargin{2};
    end
  end
    
  global HOSTNAME
  % depending on where we are running whomever calls this, the data
  % location will be different
  if isempty(HOSTNAME)
    % can set this globally (recommended if running with jvm)
    % or ask operating system (works fine when 'matlab -nojvm')
    [s,HOSTNAME] = unix('hostname');
    HOSTNAME = lower(deblank(HOSTNAME));
  end;

  reduce=0;
  switch HOSTNAME
   case {'wew.fmri.cs.cmu.edu'}
    % Wei's desktop
    dataplace = '/usr1/fMRI/Data';
    reduce = 3; % laptop specific reduced set of subjects/ROIs
   case {'human','human.learning.cs.cmu.edu'}
    % Tom's laptop
    dataplace = 'E:\fMRI\Data';
    reduce = 3; % laptop specific reduced set of subjects/ROIs    
   case {'animal','animal.learning.cs.cmu.edu'}
    % Tom's animal laptop
    dataplace = 'E:/fMRI/data';
    reduce = 3; % laptop specific reduced set of subjects/ROIs    
   case {'natural','natural.learning.cs.cmu.edu'}
    % Tom's laptop in Windows2000
    dataplace = 'F:\fMRI_data';
    reduce = 3; % laptop specific reduced set of subjects/ROIs    
   case {'pingo','pingo.pc.cs.cmu.edu'}
    % francisco's laptop
    dataplace = '/home/fpereira/projecto/FMRI/fmri_data';
    reduce = 2; % laptop specific reduced set of subjects/ROIs
   case {'liatris','liatris.cald.cs.cmu.edu'}
    % xuerui's machine
    dataplace = '/usr0/xuerui/fMRI/data';
    reduce = 3; % workstation specific reduced set of subjects/ROIs
   case {'gs2011','gs2011.sp.cs.cmu.edu'}
    % rebecca's machine
    dataplace = '/usr0/rah/FMRIdata';
    reduce = 1; % workstation specific reduced set of subjects/ROIs
    % add a comment
   case {'gs233','gs233.sp.cs.cmu.edu'}
    % stefan's machine
    dataplace = '/usr2/proiecte/fmri/fmri_data';
    reduce = 3; % workstation specific reduced set of subjects/ROIs
    % add a comment
    % hope it works :)
   case {'gs3098', 'gs3098.sp.cs.cmu.edu'}
    % indra's machine
    dataplace = '/usr0/fmri/fmri_data';
   case {'arjuna', 'arjuna.wv.cc.cmu.edu'}
    % indra's laptop
    dataplace = '/Users/irustand/fmri_data';

   case {'calla','calla.cald.cs.cmu.edu'}
    % WHOSE?? laptop in Windows2000
    dataplace = 'D:\LAB of fMRI\Raw data';
    reduce = 6; % laptop specific reduced set of subjects/ROIs
    
   case {'multivac.cald.cs.cmu.edu'}

    % MULTIVAC
    
    switch dataset
     case {'data-wordcategory','data-starplus-sp','data-starplus-ps','data-brainlex','data-brainlexNotNormalized','data-brainlexNew'}
       dataplace = '/usr1/fmri/fmri_data';
     case {'data-exemplar','data-categories','data-syntamb2','data-sixcategories','data-threecategories','data-threecategoriesNew'}
      dataplace = '/usr2/fmri/fmri_data';
    end
    
   otherwise
    
    % EVERY MACHINE ON THE CALD LAB BUT MULTIVAC
    
    switch dataset
     case {'data-wordcategroy','data-starplus-sp','data-starplus-ps','data-brainlex','data-brainlexNotNormalized','data-brainlexNew'}    
      dataplace = '/MULTIVAC/usr1/fmri/fmri_data';
     case {'data-exemplar','data-categories','data-syntamb2','data-sixcategories','data-threecategories','data-threecategoriesNew'}
      dataplace = '/MULTIVAC/usr2/fmri/fmri_data';
    end
  end
  
  % subjects available
  switch dataset

   case 'data-wordcategory'  
    subjects = {'01244B/Rob/train'};
    switch reduce
     case 1
      subjects = {'01244B/Rob/train'};
      mask     = ones(1,length(subjects));
     otherwise
      subjects = sort(subjects);
      mask     = ones(1,length(subjects));
    end       

    % ROIS not available
    rois = {'PART'};
    dimensions = [64 64 16]; nconds = 14;

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
      mask     = ones(1,length(subjects));
     case 2
      subjects = {'04847_20_sp3terC','04799_20_sp3terC'};
      mask     = ones(1,length(subjects));
     case 3
      subjects = {'04847_20_sp3terC','04820_20_sp3terC','04799_20_sp3terC'};
      mask     = ones(1,length(subjects));


     otherwise
      subjects = sort(subjects);
      mask     = ones(1,length(subjects));
      %      mask  = [0 1 0 1 1 0 0 0 0 0 1 0 0 0 1 0 1]; % trialQuality mask    
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
      mask     = ones(1,length(subjects));
     case 2
      subjects = {'04847_30_ps3terC','04799_30_ps3terC'};
      mask     = ones(1,length(subjects));
     case 3
      subjects = {'04847_30_ps3terC','04820_30_ps3terC','04799_30_ps3terC'};
      mask     = ones(1,length(subjects));
     otherwise
      subjects = sort(subjects);
      mask     = ones(1,length(subjects));
      %      mask  = [0 1 0 1 1 0 0 0 0 0 1 0 0 0 1 0 1]; % trialQuality mask    
    end

    % ROIS available
        rois = {'CALC' 'LFEF' 'LIPL' 'LIT' 'LPPREC' 'LSPL' 'LTRIA' 'RFEF' 'RIPS' 'ROPER' 'RSGA' 'RT' 'SMA' 'LDLPFC' 'LIPS' 'LOPER' 'LSGA' 'LT' 'RDLPFC' 'RIPL' 'RIT' 'RPPREC' 'RSPL' 'RTRIA'};
    
   dimensions = [64 64 8]; nconds = 3;

    
   case 'data-syntamb2'
    subjects = {'02945_22' '02953_32' '02964_22' '02999_24' '03002_31' '03003_32'};

    switch reduce
     case 1
      subjects = {'02945_22'};
      mask     = ones(1,length(subjects));
     case 4
      subjects = {'02945_22' '03003_32'};
      mask     = ones(1,length(subjects));
     otherwise
      subjects = sort(subjects);
      mask     = ones(1,length(subjects));
    end

    rois = {'LB' 'LT' 'RB' 'RT'};
    dimensions = [64 64 7]; nconds = 7;

    
   case 'data-brainlex'

    % need to fetch the detrend.dat files for other subjects
    % before extracting them
    subjects = {'08057' '08170' '08179' '08298'};    
    %subjects = {'08057' '08179' '08170' '08240' '08276' '08298' '08330'};

    switch reduce
     case 1
      subjects = {'08057'};
      mask     = ones(1,length(subjects));
     case 2
      subjects = {'08057' '08298'};
      mask     = ones(1,length(subjects));
     otherwise
      subjects = sort(subjects);
      mask     = ones(1,length(subjects));
    end
    
    % minimal set of ROIs that all the subjects have
    rois = {'ACING' 'LIES' 'LIT' 'LSES' 'LT' 'RIES' 'RIT' 'RSES' ...
	    'RT' 'SMFP' 'LDLPFC' 'LIPL' 'LOPER' 'LSGA' 'LTRIA' ...
	    'RDLPFC' 'RIPL' 'ROPER' 'RSGA' 'CALC' 'LFEF' 'LIPS' ...
	    'LPPREC' 'LSPL' 'OP' 'RFEF' 'RIPS' 'RSPL' 'SMA' 'RTRIA' ...
	    'LCBELL' 'RCBELL' 'RPPREC'};

    % There's also ALLBRAIN, with more voxels than all of these put together.
%    rois = {'ALLBRAIN'};

    dimensions = [64 64 14]; nconds = 4;
    
    
   case {'data-brainlexNotNormalized','data-brainlexNew'}

    % Both are the same as brainlex, but data was extracted from
    % binaries that were corrected for timing of acquisition of
    % slices and had missing images smoothed over (6 subjects total).
    
    % In addition, notNormalized is not normalized % above fixation.
    % For some reason, the detrend.dat binary for 08330 was not
    % correcte or smoothed, so that was not extracted.
    
    subjects = {'08057' '08170' '08179' '08240' '08276' '08298'};    
    %subjects = {'08057' '08179' '08170' '08240' '08276' '08298' '08330'};

    
    switch reduce
     case 1
      subjects = {'08057'};
      mask     = ones(1,length(subjects));
     case 2
      subjects = {'08057' '08298'};
      mask     = ones(1,length(subjects));
     otherwise
      subjects = sort(subjects);
      mask     = ones(1,length(subjects));
    end
    
    % minimal set of ROIs that all the subjects have
    rois = {'ACING' 'LIES' 'LIT' 'LSES' 'LT' 'RIES' 'RIT' 'RSES' ...
	    'RT' 'SMFP' 'LDLPFC' 'LIPL' 'LOPER' 'LSGA' 'LTRIA' ...
	    'RDLPFC' 'RIPL' 'ROPER' 'RSGA' 'CALC' 'LFEF' 'LIPS' ...
	    'LPPREC' 'LSPL' 'OP' 'RFEF' 'RIPS' 'RSPL' 'SMA' 'RTRIA' ...
	    'LCBELL' 'RCBELL' 'RPPREC'};

    % There's also ALLBRAIN, with more voxels than all of these put together.
%    rois = {'ALLBRAIN'};

    dimensions = [64 64 14]; nconds = 4;


   case 'data-categories'

    subjects = {'05483' '05502' '05549' '05675' '05688' '05695' '05797' '05813' '05833' '05838' '05868' '05874' '05886'};

    switch reduce
     case 1
      subjects = {'05886'};
      mask     = ones(1,length(subjects));
     case 2
      subjects = {'05868' '05688' '05833'};
      mask     = ones(1,length(subjects));
     otherwise
      subjects = sort(subjects);
      mask     = ones(1,length(subjects));
    end

    % some subjects don't have LCBELL, so that is out, and so is RCBELL
    %rois     = {'CALC' 'LIT' 'LT' 'LIPL' 'LIES' 'LOPER' 'LTRIA' 'LIPS' 'ACING'   'LDLPFC' 'LSGA'  'LFEF'  'LPPREC'  'LSPL' 'OP'  'SMA' 'LSES'  'SMFP'};
    %rois = {'RDLPFC'  'RIPL'  'ROPER'   'RSGA'  'RTRIA' 'RFEF'    'RIPS'  'RPPREC'  'RSPL' 'RIES'    'RIT'   'RSES'    'RT'};
    rois     = {'CALC' 'LIT' 'LT' 'LIPL' 'LIES' 'LOPER' 'LTRIA' 'LIPS' 'ACING' 'LDLPFC' 'LSGA' 'LFEF' 'LPPREC' 'LSPL' 'OP' 'SMA' 'LSES' 'SMFP' 'RDLPFC' 'RIPL' 'ROPER' 'RSGA' 'RTRIA' 'RFEF' 'RIPS' 'RPPREC' 'RSPL' 'RIES' 'RIT' 'RSES' 'RT'};
    dimensions = [64 64 16]; nconds = 7;
  

   case 'data-sixcategories'
    subjects = {'233B' '86B' '77B' '329B' '332B' '424B' '474B' '496B'};
    
    switch reduce
     otherwise
      subjects = sort(subjects);
      mask     = ones(1,length(subjects));
    end
    
    % minimal set of ROIs that all the subjects have
%rois = {'ACING' 'CALC' 'LCBELL' 'LDLPFC' 'LFEF' 'LIES' 'LIPL' 'LIPS' 'LIT' 'LOPER' 'LPPREC' 'LSES' 'LSGA' 'LSPL' 'LT' 'LTRIA' 'OP' 'RCBELL' 'RDLPFC' 'RFEF' 'RIES' 'RIPL' 'RIPS' 'RIT' 'ROPER' 'RPPREC' 'RSES' 'RSGA' 'RSPL' 'RT' 'RTRIA' 'SMA' 'SMFP'};
% removedrois = {};
        
    % There's also ALLBRAIN, with more voxels than all of these put together.
             rois = {'ALLBRAIN'};

    dimensions = [64 64 16]; nconds = 7;
    
   
   case 'data-threecategories'

    subjects = {'354B','357B','362B','367B','371B'};
    switch reduce
     otherwise
subjects = sort(subjects);
mask     = ones(1,length(subjects));
    end
    
    % minimal set of ROIs that all the subjects have
rois = {'ACING' 'CALC' 'LCBELL' 'LDLPFC' 'LFEF' 'LIES' 'LIPL' 'LIPS' 'LIT' 'LOPER' 'LPPREC' 'LSES' 'LSGA' 'LSPL' 'LT' 'LTRIA' 'OP' 'RCBELL' 'RDLPFC' 'RFEF' 'RIES' 'RIPL' 'RIPS' 'RIT' 'ROPER' 'RPPREC' 'RSES' 'RSGA' 'RSPL' 'RT' 'RTRIA' 'SMA' 'SMFP'};
    % There's also ALLBRAIN, with more voxels than all of these put together.
%        rois = {'ALLBRAIN'};

    dimensions = [64 64 16]; nconds = 4;


   case 'data-threecategoriesNew'

    %    subjects = {'354B','357B','362B','367B','371B'};
    
    %    subjects = {'354B' '362B' '371B'}; % subjects with LHG/RHG ROIs
    subjects = {'357B' '367B'}; % subjects without

    switch reduce
     otherwise
subjects = sort(subjects);
mask     = ones(1,length(subjects));
    end
    
    % minimal set of ROIs that all the subjects have
%rois = {'ACING' 'CALC' 'LCBELL' 'LDLPFC' 'LFEF' 'LIES' 'LIPL' 'LIPS' 'LIT' 'LOPER' 'LPPREC' 'LSES' 'LSGA' 'LSPL' 'LT' 'LTRIA' 'OP' 'RCBELL' 'RDLPFC' 'RFEF' 'RIES' 'RIPL' 'RIPS' 'RIT' 'ROPER' 'RPPREC' 'RSES' 'RSGA' 'RSPL' 'RT' 'RTRIA' 'SMA' 'SMFP'};
    % There's also ALLBRAIN, with more voxels than all of these put together.
    rois = {'ALLBRAIN'};
%    rois = {'LHG' 'RHG'};

    dimensions = [64 64 16]; nconds = 4;

    
  end

  information.host = HOSTNAME;
  information.dataplace = dataplace;
  information.subjects  = subjects;
  information.mask      = mask;
  information.rois      = rois;
  information.dimensions = dimensions;
  information.nconds    = nconds;
