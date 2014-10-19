%
% Computation of SVD in a way that is quicker and uses less memory
% than the svd that comes with matlab.
%


function [U,S,V] = compute_fastSVD(X)

tmp   = X * X';
[U,D] = eig(tmp);
clear tmp;
S     = sqrt(D);
tmp   = diag(S);
% [indices,order] = sort(tmp,1,'descend'); % only works in Matlab R14
[indices,order] = sort(tmp,1);
indices = flipud(indices);
order   = flipud(order);
S     = diag(indices);
U     = U(:,order);
V     = X' * U * inv(S);




function [] = testThis;

X = randn(10,100);
[U1,S1,V1] = svd(X);
minp = min(size(X,1),size(X,2));
S1 = S1(:,1:minp);
V1 = V1(:,1:minp);

[U2,S2,V2] = compute_fastSVD(X);

X1 = U1*S1*V1';
X2 = U2*S2*V2';

clf
imagesc(abs(X1-X));
colorbar
pause
clf;
imagesc(abs(X2-X));
colorbar
pause
clf;
imagesc(abs(X1-X2));
colorbar
clf;
imagesc(abs(S1-S2));
colorbar


