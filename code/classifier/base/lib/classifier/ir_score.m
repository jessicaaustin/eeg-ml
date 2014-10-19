function [precision, recall, F, DPrime] = ir_score(confusion);

%                 actual class
%                 (observation)
%
% predicted class tp               fp
% (expectation)   (correct result) (unexpected result)
%
%                 fn               tn
%                 (missing result) (correct absence of result)

if isempty(confusion) || length(confusion) == 1,
	precision = NaN;
	recall = NaN;
	F = NaN;
	DPrime = NaN;
	return;
end

tp = confusion(1,1);
fp = confusion(1,2);
fn = confusion(2,1);
tn = confusion(2,2);

precision = tp / (tp + fp);

recall = tp / (tp + fn);

F = 2 * (precision * recall) / (precision + recall);

hit = tp / (tp + fn);
falsealarm = fp / (fp + tn);
zHit = norminv(hit);
zfa = norminv(falsealarm);

DPrime = zHit - zfa;
