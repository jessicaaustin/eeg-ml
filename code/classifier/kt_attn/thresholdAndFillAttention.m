function attention = thresholdAndFillAttention(attention, attentionParams)
% given a sequence with continuous attention values, threshold 
% the attention to give "inattentive" (1) or "attentive" (2)
% for encounters with no EEG signal, fill based on the user's attention
% distribution (gaussian with mean and variance calculated based off data)

meanAttention = attentionParams.meanAttention;
varianceAttention = attentionParams.varianceAttention;

badidx = attention == -1;
attention(badidx) = varianceAttention.*randn(sum(badidx),1) + meanAttention;

attentionThresh = meanAttention;
attentiveidx = attention >= attentionThresh;
inattentiveidx = attention < attentionThresh;
attention(attentiveidx) = 2;
attention(inattentiveidx) = 1;

end
