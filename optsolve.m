function [ optsolve ] = optsolve(problem, algorithm, iterants ,kernel, blurredimage, params, endresult, conv)
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

applyD = @(x) cat(3, applyD1(x), applyD2(x));

ApplyA = @(x) cat(3, applyK(x), applyD1(x), applyD2(x));
ApplyATrans = @(x) applyKTrans(x(:,:,1)) + applyD1Trans(x(:,:,2)) + applyD2Trans(x(:,:,3)); 

applyDTrans = @(y) applyD1Trans(y(:,:,1)) + applyD2Trans(y(:, :, 2));

% Function which computes the (I + K^TK + D^TD)x where x in R^(m x n)
% matrix and the eigenvalues of I + t*t*K^TK + t*t*D^TD; here t is the
% stepsizes

% dirty, will work for now
if strcmp(algorithm, "douglasrachfordprimal")==1
    t = params.tprimaldr;
    s=t;
elseif strcmp(algorithm, "douglasrachfordprimaldual")
    t = params.tprimaldualdr;
    s=t;
elseif strcmp(algorithm, "admm")==1
    t = params.tadmm;
    s=t;
elseif strcmp(algorithm, "chambollepock")==1
    t = params.tchamb;
    s = params.schamb;
end

if (conv)
    distance=zeros(params.maxiter-1, 1);
end


applyMat = @(x) x + applyKTrans(applyK(x)) + applyDTrans(applyD(x));
eigValsMat = ones(numRows, numCols) + t*t*eigArry_KTrans.*eigArry_K + t*t*eigArry_D1Trans.*eigArry_D1...
    + s*t*eigArry_D2Trans.*eigArry_D2;

%R^(m x n) Computing (I + K^T*K + D^T*D)^(-1)*x
invertMatrix = @(x) ifft2(fft2(x)./eigValsMat); 


if strcmp(algorithm, "douglasrachfordprimal")==1
    update_iterants = @(iterants, k) PrimalDRSplit(iterants, problem, params, ApplyA, invertMatrix, ApplyATrans, (k==params.maxiter));
elseif strcmp(algorithm, "douglasrachfordprimaldual")
    update_iterants = @(iterants, k) PrimalDualDRSplit(iterants, problem, blurredimage, params, ApplyA, invertMatrix, ApplyATrans, (k==params.maxiter));
elseif strcmp(algorithm, "admm")==1
    update_iterants = @(iterants) AdmmUpdate(iterants, problem, params);
elseif strcmp(algorithm, "chambollepock")==1
    L2 = max(max(abs(eigArry_K).^2 + abs(eigArry_D1).^2 + abs(eigArry_D2).^2));
    check = params.schamb * params.tchamb * L2;
    fprintf('Stability Check (s*t*L^2): %.4f\n', check);
    update_iterants = @(iterants, k) ChambolleUpdate(iterants, problem, params, ApplyA, ApplyATrans, blurredimage);
end
tic;
for k=1:params.maxiter
    iterants = update_iterants(iterants, k);
    if mod(k, 50) == 0
        fprintf('Iter %4d | Time: %.2f s\n', k, toc);
    end
    if conv && k~=params.maxiter
        distance(k, 1)=norm((iterants.x-endresult),2);
    end
   
end
if conv
    optsolve = distance;
else 
    optsolve = iterants;
end