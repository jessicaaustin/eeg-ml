% Save the classifier
%
% This function saves the classifier which is trained by trainClassifier.
%
% In:
% - classifier (learned classifier, a structure)
% - savingPlace (the absolute path where the classifier will be saved)
% - savingFileName (the classifier will be stored in savingFileName)
%
% History
% - 06 Apr 05 - Wei Wang 
%
% Ex:
% - saveClassifier(classifier, '.','classifier.mat')

function saveClassifier( varargin )
    
    l=length(varargin);
    if l<3; help saveClassifier; return; end
    
    classifier = varargin{1};
    savePlace       = varargin{2};
    savingFileName  = varargin{3};
    
    if ~isdir(savePlace)
        fprintf([savePlace,'\nIS NOT A DIR\nPLEASE CHOOSE A NEW SAVEING DIR\n']);
        return;
    end
    
    save([savePlace,'/',savingFileName],'classifier');