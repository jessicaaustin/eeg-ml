function attention = thresholdAndFillAttention(attention, timeElapsed, params)
% given a sequence with continuous attention values, threshold 
% the attention to give "inattentive" (1) or "attentive" (2)
% for encounters with no EEG signal, fill based on the user's attention
% distribution (gaussian with mean and variance calculated based off data)

meanAttention = params.meanAttention;
varianceAttention = params.varianceAttention;

% fill in bad data
badidx = attention == -1;
attention(badidx) = varianceAttention.*randn(sum(badidx),1) + meanAttention;

timeThresh = params.meanTimeElapsed;
attentionThresh = meanAttention;

% create binary attention variable based on EEG and time elapsed
attentiveidx = (attention >= attentionThresh) & (timeElapsed < timeThresh);
attention(~attentiveidx) = 1;
attention(attentiveidx) = 2;


end
