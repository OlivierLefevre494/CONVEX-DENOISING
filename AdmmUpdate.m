function [out] = AdmmUpdate(iterants, problem, blurredimage, params, ApplyA, invertMatrix, ApplyATrans, last)
    if last
        u = iterants.u;
        y = iterants.y;
        w = iterants.w;
        z = iterants.z;

        t = params.tadmm;

        out = invertMatrix(double(u) + ApplyATrans(double(y)) - (1/t) .* (double(w) + ApplyATrans(double(z))));
    else
        
        t = params.tadmm;
        rho = params.rhoadmm;

        u = iterants.u;
        y = iterants.y;
        w = iterants.w;
        z = iterants.z;

        % x update: x_k = (I + A^T*A)^-1 (u_k-1 + A^T*y_k-1 - 1/t*(w_k-1 + A^T*z_k-1))
        x = invertMatrix(double(u) + ApplyATrans(double(y)) - (1/t) .* (double(w) + ApplyATrans(double(z))));
        Ax = ApplyA(x);

        % u update: u_k = prox_(1/t f)(p*x_k + (1 - p)*u_k + w_k/t)
        u = NormalizeImage(rho .* double(x) + (1 - rho) .* double(u) + double(w)./t);

        % y update: y_k = prox_(1/t g)(p*A*x_k + (1 - p)*y_k-1 + z_k-1/t)
        % define argument of the prox first:
        y_arg = rho .* double(Ax) + (1 - rho) .* double(y) + double(z)./t; 

        if strcmp(problem, "l2") == 1
            y = cat(3, TwoNormProx(y_arg(:,:,1), blurredimage, (1/t)), IsoProx(y_arg(:,:,2), y_arg(:,:,3), params.gammal2 / t));
        elseif strcmp(problem, "l1") == 1
            y = cat(3, OneNormProx(y_arg(:,:,1), blurredimage, (1/t)), IsoProx(y_arg(:,:,2), y_arg(:,:,3), params.gammal1 / t));
        end

        w = double(w) + t .* (double(x) - double(u));
        z = double(z) + t .* (double(Ax) - double(y));

        iterants.x = x;
        iterants.u = u;
        iterants.y = y;
        iterants.w = w;
        iterants.z = z;
        out = iterants;
    end
    iterants.output = invertMatrix(double(u) + ApplyATrans(double(y)) - (1/t) .* (double(w) + ApplyATrans(double(z))));
end
