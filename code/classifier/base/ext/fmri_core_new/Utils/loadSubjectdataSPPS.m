% loadSubjectdataSPPS(Subjects,ROIs)
%
% Loads info+data+meta for sp and ps, subjects, and ROIs from disk.
% Subjects is a cell array of two IDs of subjects.  ROIs is a cell array of
% ROI names. 
%
% Returns info data and meta, corresponding to the data loaded for corresponding subjects in sp and ps 
%
% Example:  
%  [info,data,meta] = loadSubjectdataSPPS({'04799_20_sp3terC', '04799_30_ps3terC'}, {'LIT' 'CALC'}) 
% merges starplus SP and PS datasets into a single IDM with all the
% trials in both.
%
% trial with condition       becomes a trial of condition
% SP cond 0 (noisy)          0
% SP cond 1 (fixation)       1
% SP cond 2 (affirmative)    2
% SP cond 3 (negative)       3
% PS cond 0 (noisy)          0
% PS cond 1 (fixation)       1
% PS cond 2 (affirmative)    4
% PS cond 3 (negative)       5
%
% Subjects which can be used are 
%04799, 04805, 04820, 04847, 04958, 05005, 05018, 05093, 05131,
%05675, 05680, 05695, 05710
% with total of 24 Roi's
% rois = {'CALC' 'LDLPFC' 'LFEF' 'LIPL' 'LIPS' 'LIT' 'LOPER' 'LPPREC' 'LSGA' 'LSPL' 'LT' ...
%         'LTRIA' 'RDLPFC' 'RFEF' 'RIPL' 'RIPS' 'RIT' 'ROPER' 'RPPREC' 'RSGA' 'RSPL' 'RT' 'RTRIA' 'SMA'}
% 05099
% with total of 7 Roi's
% rois = {'CALC' 'LIPL' 'LT' 'LTRIA' 'LOPER' 'LIPS' 'LDLPFC'};
% 
%
% History
% - Sep 07, 2004 - Wei (adhoc_joinStarplusDatasetsMerged from Francisco)


function [rinfo,rdata,rmeta] = loadSubjectdataSPPS( varargin )

    l = length(varargin);
  
    if l < 1
        fprintf(1,'syntax: loadSubjectdataSPPS(subjects, rois, [<use separate rois>]\n\n');
        fprintf('- subject - which subject inside the study from data-starplus-sp and data-starplus-ps\n');
        fprintf('- rois - which ROIs to use\n');
        fprintf('- useSeparateROIs (default=0) - if 1, return separate idm for each roi, else merge them all into a single idm\n');
        fprintf('\n');
        return;
    else
        subject1        = varargin{1}{1};
        subject2        = varargin{1}{2};
        rois            = varargin{2};
        useSeparateROIs = 0;
        if l > 2
            useSeparateROIs = varargin{3};
        end
        nrois           = length(rois);
    end
    
    study1 = 'data-starplus-sp';
    study2 = 'data-starplus-ps';
    
    [info1, data1, meta1] = loadSubjectdataMult(study1, subject1, rois, useSeparateROIs);
    [info2, data2, meta2] = loadSubjectdataMult(study2, subject2, rois, useSeparateROIs);
 
    [rinfos, rdatas, rmetas] = adhoc_joinStarplusDatasetsMerged(info1, data1, meta1, info2, data2, meta2);
    rinfo = rinfos{1};
    rdata = rdatas{1};
    rmeta = rmetas{1};
      
    
    
    
    
 function [rinfo,rdata,rmeta] = adhoc_joinStarplusDatasetsMerged(rinfo1,rdata1,rmeta1,rinfo2,rdata2,rmeta2)
  
    nrois = length(rmeta1);
  
    % Label all the 2/3 trials in PS as 4/5
    ntrials1 = length(rdata1{1});
    ntrials2 = length(rdata2{1});
    ntrials  = ntrials1 + ntrials2;
  
    for r=1:nrois
        for nt=1:ntrials2
            switch rinfo2{r}(nt).cond
                case {2,3}
	                if rinfo2{r}(nt).cond == 2; rinfo2{r}(nt).cond = 4; else rinfo2{r}(nt).cond = 5;end 
                otherwise
	                % leave them as they are
            end
        end
    end
  
    % concatenate them
    for r=1:nrois
        idx = 1;
        rdata{r} = cell(ntrials,1);
        
        rmeta{r} = rmeta1{r};
        rmeta{r}.ntrials = ntrials;
        rmeta{r}.study   = 'data-starplus-spps';
        
        dotDash = findstr('_', rmeta{r}.subject);
        rmeta{r}.subject   = rmeta{r}.subject(1 : dotDash(1)-1);
    
        for nt=1:ntrials1
            rinfo{r}(idx) = rinfo1{r}(nt);
            rdata{r}{idx} = rdata1{r}{nt};
            idx = idx + 1;
        end

        for nt=1:ntrials2
            rinfo{r}(idx) = rinfo2{r}(nt);
            rdata{r}{idx} = rdata2{r}{nt};
            idx = idx + 1;
        end
        
    end
    
    for r=1:nrois
        rmeta{r}.nsnapshots = sum([rinfo{r}.len]);
    end