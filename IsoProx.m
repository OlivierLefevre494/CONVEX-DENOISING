function [ out ] = IsoProx(x1, x2, t)
    
    % Compute the magnitude of each pixel
    magnitude = sqrt(x1.^2 + x2.^2);
    
    % zero if magnitude <= t, otherwise shrink toward zero
    alpha = max(1 - t./magnitude, 0);
    alpha(magnitude == 0) = 0;  % Avoid division by zero
    
    % Return the shrunk gradient stacked as 3D array
    out = cat(3, alpha .* x1, alpha .* x2);
end