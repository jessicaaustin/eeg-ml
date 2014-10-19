function Y = cellsplit(X)

Y = {};
for i = 1:length(X),
	Y{i} = X(i); 
end
