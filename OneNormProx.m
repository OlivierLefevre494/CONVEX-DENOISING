function [ out ] = OneNormProx(x1, b, t)
    % 1. Calculate the difference
    diff = x1 - b;
    
    % sign(diff) preserves the direction
    % max(abs(diff) - t, 0) shrinks the magnitude
    res = sign(diff) .* max(abs(diff) - t, 0);
    
    % Add the background back to get the final result
    out = res + b;
end