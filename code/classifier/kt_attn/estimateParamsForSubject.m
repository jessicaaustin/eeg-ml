function [LL,p,A,B,C,D, ...
    l0, p_learn, p_forget, p_guess, p_slip, p_learn_a, p_dontlearn_a, p_guess_a, p_slip_a] ...
    = estimateParamsForSubject(sequences, Oidx, model)

%% EM params

max_iter = 500;
thresh = 1e-8;
verbose = 1;

%% Initial estimates for KT-Attn params

% a complete guess for the learning rate, and assume that the forgetting
% rate is near zero
p0 = [0.5; 
      0.5];
p_learn_0 = 0.7;
p_forget_0 = 0.01;
A0 =[1-p_forget_0  p_forget_0
     1-p_learn_0   p_learn_0];
p_guess_0 = 0.4;
p_slip_0 = 0.2;
B0 = [p_guess_0   1-p_guess_0
      p_slip_0    1-p_slip_0];
p_learn_a_0 = 0.7;
p_dontlearn_a_0 = 0.8;
C0 = [p_dontlearn_a_0   1-p_dontlearn_a_0
      1-p_learn_a_0     p_learn_a_0];
p_guess_a_0 = 0.6;
p_slip_a_0 = 0.7;
D0 = [p_slip_a_0      1-p_slip_a_0
      1-p_guess_a_0   p_guess_a_0];
       
  
%% Put observations into cell array

O = length(Oidx);
x = {};

for o=Oidx(:)'
    
    ASRobservations = sequences.accept{o};
    
    if strcmp(model, 'KTAttn')
        EEGobservations = thresholdAndFillAttention(sequences.attention{o}, sequences);
        x{end+1} = [ASRobservations';
            EEGobservations'];
    elseif strcmp(model, 'KT')
        x{end+1} = ASRobservations';
    else
        error('model %s not found!\n', model);
    end
    
end

%% E-M algorithm

if strcmp(model, 'KTAttn')
    [LL, p, A, B, C, D] = learn_dhmm_ktattn(x, p0, A0, B0, C0, D0, max_iter, thresh, verbose);
elseif strcmp(model, 'KT')
    [LL, p, A, B] = learn_dhmm(x, p0, A0, B0, max_iter, thresh, verbose);
    C = ones(2,2);
    D = ones(2,2);
else
    error('model %s not found!\n', model);
end

%% pull out model params

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