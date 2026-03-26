function [ out ] = IsoProx(x1,x2,t)
% Computes the proximal operator of the IsoNorm at cat (x1, x2) with t = t,
% x1 and x2 should both! have shape R^(nxn)
[numRows, numCols] = size(x1);
res = cat(3, x1, x2);

for i = 1:numRows
    for j = 1:numCols
        if ((x1(i,j,1)^2 + x2(i,j,1)^2)^0.5 > t)
            alpha = 1.0-(t/(x1(i,j,1)^2 + x2(i,j,1)^2)^0.5);
        else 
            alpha = 0.0;
        end
        res(i, j, 1) = alpha*x1(i,j,1);
        res(i,j, 2) = alpha*x2(i,j,1);
    end
end

out = res(:,:,1:2);
end