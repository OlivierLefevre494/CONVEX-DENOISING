function [ optsolve ] = optsolve(problem, algorithm, iterants ,kernel, blurredimage, params, endresult, conv)
arguments
    problem % l1 or l2 problem
    algorithm % name of the algorithm
    iterants % initialization of the iterants for the algorithm (can use DefaultInitializeIterants.m)
    kernel % blurring kernel
    blurredimage % corrupted image as an mxn matrix
    params % algorithm params (can use DefaultParams.m)
    endresult = 0 % expected end result for algorithm (if evaluating convergence)
    conv = false % convergence study mode (on/off)

    % By default, convergence study mode is not activate.
end

[numRows, numCols] = size(blurredimage); %numRows = m, numCols = n

%computes the numRow x numCol matrix of the eigenvalues for K and D1 and
%D2; Here D1 = I oplus D1 in the paper and D2 = D1 oplus I.
eigArry_K = eigValsForPeriodicConvOp(kernel, numRows, numCols);
eigArry_D1 = eigValsForPeriodicConvOp([-1,1]', numRows, numCols);
eigArry_D2 = eigValsForPeriodicConvOp([-1,1], numRows, numCols);

%computes numRow x numCol matrix of the eigenvalues for K^T and D1^T and
%D2^T;
eigArry_KTrans = conj(eigArry_K);
eigArry_D1Trans = conj(eigArry_D1);
eigArry_D2Trans = conj(eigArry_D2);

%Functions which compute Kx, D1x, D2x, Dxt, K^Tx, D1^Tx, D2^Tx, and D^Ty.
%Note for all the x functions, the input x is in R^(m x n) and outputs into
%R^(m x n) except for D which outputs into 2 concat. R^(m x n) matrices;
%For D^Ty, y is two m x n matrices concatanated and outputs into R^(m x n)
applyK = @(x) applyPeriodicConv2D(x, eigArry_K);
applyD1 = @(x) applyPeriodicConv2D(x, eigArry_D1);
applyD2 = @(x) applyPeriodicConv2D(x, eigArry_D2);

applyKTrans = @(x) applyPeriodicConv2D(x, eigArry_KTrans);
applyD1Trans = @(x) applyPeriodicConv2D(x, eigArry_D1Trans);
applyD2Trans = @(x) applyPeriodicConv2D(x, eigArry_D2Trans);

ApplyA = @(x) cat(3, applyK(x), applyD1(x), applyD2(x));
ApplyATrans = @(x) applyKTrans(x(:,:,1)) + applyD1Trans(x(:,:,2)) + applyD2Trans(x(:,:,3)); 

switch algorithm
    case 'douglasrachfordprimal'
        t = params.tprimaldr;
    case 'douglasrachfordprimaldual'
        t = params.tprimaldualdr;
    case 'admm'
        t = params.tadmm;
    case 'chambollepock'
        t = params.tchamb;
        s = params.schamb;
    otherwise
        error('Not valid algorithm "%s"', algorithm);
end

% Resolvent for Algo 1&3
eigA2 = abs(eigArry_K).^2 + abs(eigArry_D1).^2 + abs(eigArry_D2).^2;

eigVals_noStep = 1 + eigA2;
invertMatrixNoStep = @(x) ifft2(fft2(x) ./ eigVals_noStep);
 
% Resolvent for Algo 2
eigVals_withStep = 1 + t*t*eigA2;
invertMatrixWithStep = @(x) ifft2(fft2(x) ./ eigVals_withStep);

switch algorithm
    case 'douglasrachfordprimal'
        update_iterants = @(it, k) PrimalDRUpdate(it, problem, blurredimage, params, ApplyA, invertMatrixNoStep, ApplyATrans, (k == params.maxiter));
 
    case 'douglasrachfordprimaldual'
        update_iterants = @(it, k) PrimalDualDRUpdate(it, problem, blurredimage, params, ApplyA, invertMatrixWithStep, ApplyATrans, (k == params.maxiter));
 
    case 'admm'
        update_iterants = @(it, k) AdmmUpdate(it, problem, blurredimage, params, ApplyA, invertMatrixNoStep, ApplyATrans, (k == params.maxiter));
 
    case 'chambollepock'
        L2 = max(eigA2(:));
        check = s * t * L2;
        fprintf('Stability check (s*t*||A||^2): %.4f   (must be < 1)\n', check);
        update_iterants = @(it, k) ChambolleUpdate(it, problem, blurredimage, params, ApplyA, ApplyATrans);
end

if (conv)
    distance=zeros(params.maxiter-1, 1);
end

tic;
for k=1:params.maxiter
    iterants = update_iterants(iterants, k);
    if mod(k, 50) == 0
        fprintf('Iter %4d | Time: %.2f s\n', k, toc);
    end
    if conv && k~=params.maxiter
        distance(k, 1) = norm(iterants.x(:) - endresult(:), 2);
    end
   
end
if conv
    optsolve = distance;
else 
    optsolve = iterants;
end
