% Print out a list of studies, and where data will be sought for each on
% current host machine
%
% Example: (note it takes no arguments)
%  locateStudies
%
% WISH LIST: make it look in the directories and display which subjects
% and ROIs have data on this machine

function  locateStudies()
  studies = {'data-brainlex' 'data-brainlexNew' 'data-brainlexNotNormalized' ...
      'data-categories' 'data-sixcategories'...
      'data-starplus-ps' 'data-starplus-sp' 'data-starplus-spps' 'data-syntamb2'...
      'data-threecategories'   'data-fourcategories' 'data-biEng' 'data-biPort'};
      %'data-brainlex-org' 'data-exemplar' 'data-goldcategories'...
      %'data-threecategories-old' 'data-twocategories'};

  for s = 1:length(studies)
    fprintf('%s\n',studies{s});
    information = localInformation(studies{s});   
    fprintf('\t%s/%s\n',information.dataplace,studies{s});
  end
