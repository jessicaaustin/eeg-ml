function [n, accuracy, p] = evaluate_results(gold, results)
    [acc, racc, eY, rank, confusion] = classifier_score(gold, results);
    confusion
    correct = length(find(acc));
    incorrect = length(~find(acc));
    total = length(eY);
    chi2 = chi_squared_sig_test([correct, incorrect],   [total * 0.5, total * 0.5], 0.05);
    sig_indicator = {'', '*'};

	n = length(acc);
	accuracy = length(find(acc)) / length(acc);
	p = chi2.P;

    fprintf('n=%d accuracy=%1.2f%s p=%1.2f\n', n, accuracy, sig_indicator{chi2.H+1}, p);
end
