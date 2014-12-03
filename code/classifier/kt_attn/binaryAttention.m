function attention = binaryAttention(attention)
% given a sequence with continuous attention values, threshold 
% the attention to give "inattentive" (1) or "attentive" (2)

load('attentionThresh.mat')

badidx = attention == -1;
attentiveidx = attention >= attentionThresh;
inattentiveidx = attention < attentionThresh;
attention(attentiveidx) = 2;
attention(inattentiveidx) = 1;
% attention(badidx) = -1;  % TODO: add this back in!!!!!!!

end