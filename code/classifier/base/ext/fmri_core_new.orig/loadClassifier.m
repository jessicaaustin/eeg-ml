% load the saved classifier
%
% This function loads the classifier which was trained before.
%
% In:
% - savingPlace (the absolute path where the classifier was saved)
% - savingFileName (the classifier was stored in savingFileName)
%
% Out:
% - outputClassifier
%
% History
% - 06 Apr 05 - Wei Wang 
%
% Ex:
% - [classifier]=loadClassifier('.','classifier.mat');

function [outClassifier]=loadClassifier( varargin )
    
    l=length(varargin);
    if l<2; help loadClassifier; return; end
    
    savePlace       = varargin{1};
    savingFileName  = varargin{2};
    
    if ~isdir(savePlace)
        fprintf([savePlace,'\nIS NOT A DIR\nPLEASE CHECK SAVEING DIR\n']);
        classifier=[];
        return;
    end
    
    load([savePlace,'/',savingFileName]);
    outClassifier = classifier;
    