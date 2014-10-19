% plotSMLRweights(meta,classifier_SMLR,classNum)
%
% Create a image of the classNum's learned weights which are learned by
% SMLR. minActiv and maxActiv determine the minimum and
% maximum activation levels displayed in the color scale of the image.
% minActiv and maxActiv were detemined by the min and max values of the
% weights over the weights.
%
% images are laid out in format: 
% 1  2  3  4 
% 5  6  7  8 
% 9 10 11 12
%13 14 15 16
%
% Example: display third class of weights.
%  img = plotSMLRweights(meta,classifier,3);  
% 
% History
% - 06/15/2005 Wei - created.

function [img]=plotSMLRweights(varargin)

l = length(varargin);
  
  if l < 3
    fprintf('syntax: plotSMLRweights(meta, classifier, classNum)\n');
    return;
  end
  meta = varargin{1};
  classifier = varargin{2};
  classNum = varargin{3};
  
  [weightInfo, weightData, weightMeta, minWeight, maxWeight]= transformIDM_SMLR_weightsToIDM(meta, classifier);
  img=plotSnapshot3D(weightInfo, weightData, weightMeta,classNum,1,minWeight, maxWeight);


function [rinfo, rdata, rmeta, minData, maxData]= transformIDM_SMLR_weightsToIDM(varargin)
    
l = length(varargin);
  
  if l < 2
    fprintf('syntax: transformIDM_SMLR_weightsToIDM(meta, classifier)\n');
    return;
  end
    
  meta = varargin{1};
  classifier = varargin{2};
  
  models = classifier.models;
  lastpos = size(models,1);
  weights_all=models{lastpos-1}.uncombinedWeights;
  weights=weights_all(2:end,:);
  
  minData=min(min(weights));
  maxData=max(max(weights));
  
  [numRows, numCols]=size(weights);
  
  rdata=cell(numCols,1);
  for i = 1:length(rdata)
      rdata{i}=weights(:,i)';
      rinfo(i).cond=i;
      rinfo(i).len=1;
  end
  
  rmeta=meta;
  rmeta.nsnapshots=numCols;
  rmeta.ntrials=numCols;
  