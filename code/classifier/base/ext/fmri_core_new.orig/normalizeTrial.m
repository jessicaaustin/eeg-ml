% Returns array containing normalized values of the voxels in trial Values
% are normalized in two possible ways:
%
% 1. 'perVoxel' so that the mean activity of each voxel is adjusted to zero
%
% 2. 'perBrain' so that the relative values of voxels are kept the same, but
% the mean activity over all voxels and snapshots equals zero.
%
% There is no normalization of the variance

function newTrial = normalizeTrial(trial)
  
  normalizeMethod='perVoxel'; % should be an optional input argument

  newTrial=zeros(size(trial));
  
  % calculate row vector of means, one for each voxel
  means = mean(trial);
  
  switch normalizeMethod
   case 'perVoxel'
    % normalize to set each voxel mean to zero during this trial
    for j=1:1:size(trial,1)  % for each row
      newTrial(j,:)=trial(j,:)-means;
    end

   case 'perBrain'
    % add same offset to each voxel, to set mean over all voxels & trials to
    % zero
    meanMean=mean(means)
    newTrial = trial - meanMean;
  end
  
