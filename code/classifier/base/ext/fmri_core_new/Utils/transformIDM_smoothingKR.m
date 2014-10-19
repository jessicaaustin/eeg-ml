% Smooths the IDM temporally by performing Kernel Regression with
% gaussian kernels, on a trial by trial basis
%
% In:
% - info+data+meta
% - optionally:
%   - conditions to smooth
%   - Smoothing bandwidth to use in each voxel (a 1 x #voxels vector);
%     the bandwidth is the standard deviation of the gaussian kernel.
%     By default the code will perform leave-1-point out
%     cross-validation to determine the optimal bandwidth
%
%     WARNING: for many voxels, the residuals after smoothing
%     (data - smoothed data) are strongly correlated (i.e. they
%     do not look like a sample from a mean 0 gaussian). This leads
%     to a few voxels being undersmoothed, as fitting the noise in
%     points around the one being left out will lead to a good fit
%     there as well. It's not a major issue, in that most voxels
%     seem appropriately smooth, and the worse that can happen
%     is that a voxel will be undersmoothed, rather than
%     oversmoothed (which could destroy information).
%
%     I'm in touch with someone who has written about ways of
%     correcting this, but for the time being I have no other way
%     of dealing with it.
%     
%
% Out:
% - info+data+meta, smoothed
%



% History:
% - 2004 10 06 - fpereira - created
% - 2004 10 11 - fpereira - corrected small bug in argument processing
% - 2005 05 07 - indra - add statement to clear predAcrossHCV
% Notes:
%
% See WARNING above.

function [info,data,meta] = transformIDM_smoothingKR( varargin )

%
% Process arguments
%

l = length(varargin);
if l < 3 
  help transformIDM_smoothingKR;
  return;
else
  info = varargin{1};
  data = varargin{2};
  meta = varargin{3};
  
  trialSequence = [info(:).cond];
  trialLens     = [info(:).len];
  taskTrials    = find(trialSequence > 1);
  conditions    = unique(trialSequence(taskTrials));
  nTrials       = length(trialSequence);
  nVoxels       = size(data{1},2);
  
  kernel    = [0.1 0.2 0.4 0.2 0.1];
  kernelLen = length(kernel);
  
  % possibly override list of conditions to smooth
  if l > 3; conditions = varargin{4}; end
  condMask = zeros(1,length(unique(trialSequence(taskTrials)))+1);
  condMask(conditions) = 1;
  
  % possibly fix bandwidth for each voxel
  estimateBandwidth = 1; 
  if l > 4; hVoxel = varargin{5}; estimateBandwidth = 0; end
      
  % test a few things
  taskTrialLens = trialLens(taskTrials);
  if sum(taskTrialLens < kernelLen)
    fprintf('error: all trials must be at least as long as the kernel\n');
    return;
  end
end
  
fprintf('transformIDM_smoothingKR: smoothing trials\n\t');

for nt=1:nTrials
  
  cond = info(nt).cond;
  len  = info(nt).len;
  
  % should this condition be smoothed
  if cond & condMask(cond)
    fprintf('%d ',cond);
    
    % optionally estimate the ideal bandwidth for each voxel
    
    if estimateBandwidth

      rangeOfH     = 1:0.25:3;
      predAcrossH  = zeros(len,nVoxels,length(rangeOfH));
      distAcrossH  = zeros(len,nVoxels,length(rangeOfH));
      scoresH      = zeros(length(rangeOfH),nVoxels);
%      scoresHL1O   = zeros(length(rangeOfH),nVoxels);
      totalAcrossH = zeros(len,nVoxels,length(rangeOfH));
      clear predAcrossHCV;
      
      % compute the J(h) estimate for each value of h,
      % at all voxels simultaneously
      for k = 1:length(rangeOfH)
	
	h      = rangeOfH(k);
	hVoxel = repmat(h,1,nVoxels);
	hToUse = repmat(hVoxel,len,1);
	
	% compute the predicted r(x), and also the other
	% term required by the shortcut formula
	xvalues = repmat( (1:len)', [1,nVoxels]);

	for t = 1:len
	  % current value of x, x_t
	  weights = repmat( t, [len,nVoxels]);
	  	  
	  % K( (x_t-x_i)/h ), for all i at once
	  weights = (weights - xvalues);
	  weights = weights ./ hToUse;
	  weights = normpdf(weights);
%	  weightsL1O = weights;
	  weightsCV  = weights;
	  
	  % sum of that over all i
	  total   = sum(weights,1);
	  % K( (x_t - x_i)/h ) / S_j K( (x_t - x_j)/h )
	  weights = weights ./ repmat(total, [len,1]);
	  
	  % r(x_t)
	  predAcrossH(t,:,k) = sum(data{nt} .* weights,1);

	  if 1
	  % r(x_t) leave 1 out prediction - leaves out x_t,y_t
	  % knock out one of the weights 
	  
	  %	  weightsL1O(t,:) = [];
	  %	  totalL1O   = sum(weightsL1O,1);
	  %	  weightsL1O = weightsL1O ./  repmat(totalL1O, [len-1,1]);
	  %	  tmp = data{nt};
	  %	  tmp(t,:) = [];	  
	  %	  predAcrossHL1O(t,:,k) = sum(tmp .* weightsL1O,1);

	  % r(x_t) leave a few out prediction - leaves out x_t,y_t +/-window
	  window = [-1 0 1];
	  %	  window = [-2 -1 0];
          window = window + t;
	  window = min(window,[+Inf,t,len]); % max window is len
	  window = max(window,[1,t,-Inf]); % min window is 1
%	  window = max(window,[1,1,t]);
	  if rem(t,2); window = 2:2:len; else; window = 1:2:len; end
	  
	  weightsCV(window,:) = [];
	  totalCV   = sum(weightsCV,1);
	  weightsCV = weightsCV ./  repmat(totalCV, [len-length(unique(window)),1]);
	  tmp = data{nt};
	  tmp(window,:) = [];
	  predAcrossHCV(t,:,k) = sum(tmp .* weightsCV,1);
	  end
	  
	  % the K(x_i - x_j) term on 5.20
	  distAcrossH(t,:,k) = (1 - (normpdf(0)./total)).^(-2);
	
	  totalAcrossH(t,:,k) = total;
	  
	end
	
	term = ((data{nt}-predAcrossH(:,:,k)).^2) .* distAcrossH(:,:,k);
	scoresH(k,:) = sum(term,1);
	%	scoresHL1O(k,:) = sum((data{nt}-predAcrossHL1O(:,:,k)).^2,1);
	scoresHCV(k,:) = sum( (data{nt}-predAcrossHCV(:,:,k)).^2,1);
      end
      
%      [minH,minPos] = min(scoresHCV,[],1);
      [minH,minPos] = min(scoresH,[],1);
      hVoxel = rangeOfH(minPos);
    else
      % just use the hVoxel provided
    end
    
    newData = zeros(size(data{nt}));

    for t = 1:len
      mask    = repmat( (1:len)', [1,nVoxels]);
      weights = repmat( t, [len,nVoxels]);
      weights = (weights - mask);
      hToUse  = repmat(hVoxel,len,1);
      weights = weights ./ hToUse;
      weights = normpdf(weights);
      total   = sum(weights,1);
      weights = weights ./ repmat(total, [len,1]);
      
      newData(t,:) = sum(data{nt} .* weights,1);
      if 0
	% create a convolution matrix
	row  = zeros(1,len);
	row(1:kernelLen) = kernel;
	tmp = triu(toeplitz(row));
	convMatrix = tmp(1:len-kernelLen+1,:);
	
	% convolve the trial with it
	data{nt} = convMatrix * data{nt};    
	info(nt).len = len-kernelLen+1;
      end
    end

    %
    % PLOTS FOR DEBUGGING
    %
    
    
    if 0
      % plot the voxels with extreme values of h
      hist(hVoxel)
      pause
      indices = find(hVoxel == 0.5);
      for l = 1:length(indices)
	clf;
	hold on;
	plot(data{nt}(:,indices(l)),'b');
	plot(newData(:,indices(l)),'r');
	hold off;
	pause
      end
    end

    if 0
      % plot the h curves for all voxels
      nrows = ceil(sqrt(nVoxels));
      ncols = nrows;
      for v = 1:nVoxels
	subplot(nrows,ncols,v);
	hold on;
	plot(scoresH(:,v),'b');
%	plot(scoresHL1O(:,v),'r');
	set(gca,'XTick',[]);set(gca,'YTick',[]);
	hold off;
      end
      pause
    end

    if 0
      % plot the raw/smooth curves for all voxels
      nrows = ceil(sqrt(nVoxels));
      ncols = nrows;
      for v = 1:nVoxels
	subplot(nrows,ncols,v);
	hold on;
	plot(data{nt}(:,v),'b');
	plot(newData(:,v), 'r');
	axis([0 Inf -3 3]);
	set(gca,'XTick',[]);set(gca,'YTick',[]);
	hold off;
      end
      pause
      clf;
    end
    
    if 0
      % plot the raw/smooth curves with background coloured
      % according to test results
      nrows = ceil(sqrt(nVoxels));
      ncols = nrows;
      residuals = data{nt} - newData;
      autocorr  = computeAutoCorrelation( residuals );
      colours = colormap('bone');
      
      for v = 1:nVoxels
	subplot(nrows,ncols,v);
	hold on;

	[h,p]  = kstest(residuals(:,v));
%	if p > 0.05; p = 1; end
	pos    = floor(p*64)+1;
	colour = colours(pos,:);
	
	plot(data{nt}(:,v),'b');
	plot(newData(:,v), 'r');
	axis([0 Inf -3 3]);
	set(gca,'Color',colour);
	
	set(gca,'XTick',[]);set(gca,'YTick',[]);
	hold off;
      end
      subplot(nrows,ncols,v+1);
      imagesc(colours(:,1));
      
      pause
      clf;
    end
    
    if 0
      % plot the temporal autocorrelation of residuals
      nrows = ceil(sqrt(nVoxels));
      ncols = nrows;
      residuals = data{nt} - newData;
      autocorr  = computeAutoCorrelation( residuals );
      
      for v = 1:nVoxels
	subplot(nrows,ncols,v);
	hold on;

	plot(autocorr(1:5,v));
	plot(0:5,0,'k.');
	axis([0 Inf -1 1]);
	
	set(gca,'XTick',[]);set(gca,'YTick',[]);
	hold off;
      end
      pause
      clf;
    end

    
    if 0
      % relate voxel activation level to h score
      mv = mean(data{nt},1);
      [dummy,order] = sort(mv);
       plot(dummy,hVoxel(order),'b.');
       pause
    end
    
    data{nt} = newData;
  end
  

end
fprintf('\ntransformIDM_smoothingKR: done\n');
