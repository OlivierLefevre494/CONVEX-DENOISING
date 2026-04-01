function [ out ] = TwoNormProx(x1, b, t)
% Computes the proximal operator of the 2 norm  squared at x1
[numRows, numCols] = size(x1);

out = (1/((1/t)-2))*((1/t)*(x1+b)-2*b);
end