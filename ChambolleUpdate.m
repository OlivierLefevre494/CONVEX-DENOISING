function iterants = ChambolleUpdate(iterants, problem, params, ApplyA, ApplyATrans, b)
    % unpack the iterants
    x = iterants.x;
    y = iterants.y; 
    z = iterants.z;
    s = params.schamb;
    t = params.tchamb;
    
    % gamma based on the problem type
    if strcmp(problem, 'l1')
        gamma = params.gammal1;
    else
        gamma = params.gammal2;
    end
    
    % Apply A to z
    Az = ApplyA(z); 
    
    % uphill step
    y_temp = y + s * Az;
    
    % Proximal for (y(:,:,1)) using Moreau's identity
    if strcmp(problem, 'l1')
        y(:,:,1) = y_temp(:,:,1) - s * (OneNormProx((y_temp(:,:,1) / s), b, 1/s));
    elseif strcmp(problem, 'l2')
        y(:,:,1) = y_temp(:,:,1) - s * (TwoNormProx((y_temp(:,:,1) / s), b, 1/s));
    end
    
    % 1. Evaluate the proximal operator of g/s
    prox_out = IsoProx(y_temp(:,:,2) / s, y_temp(:,:,3) / s, gamma / s);
    
    % Apply Moreau's identity
    y(:,:,2:3) = y_temp(:,:,2:3) - s * prox_out;
    
    % Apply A transpose to y
    ATy = ApplyATrans(y);
    % downhill step
    x_temp = x - t * ATy;
    % Proximal of indicator function f
    x_new = min(max(real(x_temp), 0), 1);

    z_new = 2 * x_new - x;

    % ready the iterants for the next iteration
    iterants.x = x_new;
    iterants.y = y;
    iterants.z = z_new;
end