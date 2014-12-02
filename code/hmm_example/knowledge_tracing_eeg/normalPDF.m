function [val] = normalPDF(mean, variance, x)
    val = 1./sqrt(variance*2*pi) .* exp( -1*(x-mean).^2 ./ (2*variance) );
end