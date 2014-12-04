function [LL,p,A,B,C,D, ...
    l0, p_learn, p_forget, p_guess, p_slip, p_learn_a, p_dontlearn_a, p_guess_a, p_slip_a] ...
    = estimateParamsForSubject(sequences,p0,A0,B0,C0,D0)

% EM params
max_iter = 500;
thresh = 1e-8;
verbose = 1;

% put observations into cell array
O = length(sequences.accept);
x = cell(O,1);
for o=1:O
    ASRobservations = sequences.accept{o};
    EEGobservations = thresholdAndFillAttention(sequences.attention{o}, sequences);
    x{o} = [ASRobservations';
        EEGobservations'];
end

% E-M algorithm
[LL, p, A, B, C, D] = learn_dhmm_ktattn(x, p0, A0, B0, C0, D0, max_iter, thresh, verbose);

% pull out KT-Attn params
l0 = p(2);
p_learn =  A(2,2);
p_forget = A(1,2);
p_guess = B(1,1);
p_slip = B(2,1);
p_learn_a = C(2,2);
p_dontlearn_a = C(1,1);
p_guess_a = D(2,2);
p_slip_a = D(1,1);


end