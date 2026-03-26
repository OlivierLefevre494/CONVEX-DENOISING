function [ out ] = TwoNormProx(x1, b, t)
% Computes the proximal operator of the 2 norm  squared at x1
[numRows, numCols] = size(x1);
res = zeros(numRows, numCols);

for i = 1:numRows
    for j = 1:numCols
        if (x1(i, j)-b(i,j) > t)
            res(i,j)=x1(i,j)-b(i,j)-t;
        elseif (-t<x1(i,j)-b(i,j)) && (x1(i,j)-b(i,j)<t)
            res(i,j) = 0;
        else
            res(i,j) = x1(i,j) - b(i,j)+ t;
        end

    end
end

out = res+double(b);
end