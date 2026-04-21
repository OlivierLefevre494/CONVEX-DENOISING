clear; clc; close all;
rng(42);

outdir = 'results_hyperparam';
if ~exist(outdir,'dir'); mkdir(outdir); end

% Load Image
I = im2double(rgb2gray(imread('testimages/cameraman.jpg')));
I = imresize(I, [128, 128]);
[M, N] = size(I);

padPx        = 10;
SEARCH_ITERS = 500;
FINAL_ITERS  = 1000;

algos = {'douglasrachfordprimal','douglasrachfordprimaldual','admm','chambollepock'};

kernels = { fspecial('gaussian',[9,9],1.5), fspecial('motion',9,0) };
kernelLabels = {'gauss9_s15','motion9'};

noiseL1 = @(x) imnoise(x,'salt & pepper',0.05);
noiseL2 = @(x) imnoise(x,'gaussian',0,0.01);

% Grid Search Area
tGrid     = [0.1, 0.3, 0.7, 1.0, 3.0, 8.0, 16.0];
rhoGrid   = [0.1, 0.25, 0.5, 0.8, 1.0, 1.5, 3];
sGridCP   = [0.05, 0.15, 0.35];
tGridCP   = [0.05, 0.15, 0.35];
gammaGrid = [0.001, 0.005, 0.01, 0.05, 0.1];

DEFAULT_GAMMA = 0.05;
problems = {'l1','l2'};

% Pre-build blurred+padded images
datasets = struct();
for pi_ = 1:numel(problems)
    problem = problems{pi_};
    if strcmp(problem,'l1'), nfn = noiseL1; else, nfn = noiseL2; end
    for ci = 1:numel(kernels)
        k = kernels{ci} / sum(kernels{ci}(:));
        eK = eigValsForPeriodicConvOp(k, M, N);
        b  = real(ifft2(fft2(I) .* eK));
        b  = min(max(nfn(b), 0), 1);
        b_pad = padarray(b, [padPx padPx], 'replicate', 'both');
        crop  = [padPx+1, padPx+M, padPx+1, padPx+N];
        datasets.(problem).(kernelLabels{ci}) = struct( ...
            'b_pad', b_pad, 'crop', crop, 'kernel', k, 'b', b);
    end
end

% Search loop
results = struct('problem',{},'case',{},'algo',{}, ...
                 't',{},'rho_or_s',{},'gamma',{}, ...
                 'psnr_search',{},'psnr_final',{},'ssim_final',{});

for pi_ = 1:numel(problems)
    problem = problems{pi_};
    for ci = 1:numel(kernels)
        caseLabel = kernelLabels{ci};
        ds = datasets.(problem).(caseLabel);

        for ai = 1:numel(algos)
            algo = algos{ai};
            fprintf('Tuning %s / %s / %s ', problem, caseLabel, algo);

            % Stage 1: sweep step sizes at default gamma
            [bestT, bestRoS] = stage1(algo, problem, ds, ...
                tGrid, rhoGrid, sGridCP, tGridCP, DEFAULT_GAMMA, SEARCH_ITERS);

            % Stage 2: sweep gamma at best step sizes
            [bestG, bestSearchPSNR] = stage2(algo, problem, ds, ...
                gammaGrid, bestT, bestRoS, SEARCH_ITERS);

            % Final run at FINAL_ITERS
            [finalP, finalS] = oneRun(algo, problem, ds, ...
                bestG, bestT, bestRoS, FINAL_ITERS);

            fprintf(' done  (search=%.2f  final=%.2f)\n', bestSearchPSNR, finalP);

            results(end+1) = struct( ...
                'problem', problem, 'case', caseLabel, 'algo', algo, ...
                't', bestT, 'rho_or_s', bestRoS, 'gamma', bestG, ...
                'psnr_search', bestSearchPSNR, ...
                'psnr_final',  finalP, 'ssim_final', finalS); %#ok<SAGROW>
        end
    end
end

% Summary table
fprintf('\nSUMMARY\n');
fprintf('%-4s  %-11s  %-26s  %7s  %7s  %7s  %7s  %7s\n', ...
        'obj','case','algo','t','rho/s','gamma','PSNR','SSIM');
fprintf('%s\n', repmat('-',1,88));
for k = 1:numel(results)
    r = results(k);
    fprintf('%-4s  %-11s  %-26s  %7.3g  %7.3g  %7.4g  %7.3f  %7.4f\n', ...
        r.problem, r.case, r.algo, r.t, r.rho_or_s, r.gamma, r.psnr_final, r.ssim_final);
end
fprintf('%s\n', repmat('-',1,88));

% helpers

function [bestT, bestRoS] = stage1(algo, problem, ds, ...
        tGrid, rhoGrid, sGridCP, tGridCP, gamma, iters)
    bestT = NaN; bestRoS = NaN; bestP = -inf;
    if strcmp(algo,'chambollepock')
        for si = 1:numel(sGridCP)
            for ti = 1:numel(tGridCP)
                p = oneRun(algo, problem, ds, gamma, tGridCP(ti), sGridCP(si), iters);
                fprintf('.');
                if p > bestP, bestP = p; bestT = tGridCP(ti); bestRoS = sGridCP(si); end
            end
        end
    else
        for ti = 1:numel(tGrid)
            for ri = 1:numel(rhoGrid)
                p = oneRun(algo, problem, ds, gamma, tGrid(ti), rhoGrid(ri), iters);
                fprintf('.');
                if p > bestP, bestP = p; bestT = tGrid(ti); bestRoS = rhoGrid(ri); end
            end
        end
    end
end

function [bestG, bestP] = stage2(algo, problem, ds, ...
        gammaGrid, bestT, bestRoS, iters)
    bestG = NaN; bestP = -inf;
    for gi = 1:numel(gammaGrid)
        p = oneRun(algo, problem, ds, gammaGrid(gi), bestT, bestRoS, iters);
        fprintf('.');
        if p > bestP, bestP = p; bestG = gammaGrid(gi); end
    end
end

function [ps, ss] = oneRun(algo, problem, ds, gamma, t, rhoOrS, iters)
    ps = -inf; ss = -inf;
    try
        params = DefaultParams();
        params.maxiter = iters;
        params = setGamma(params, problem, gamma);
        params = setT(params, algo, t);
        if strcmp(algo,'chambollepock')
            params.schamb = rhoOrS;
        else
            params = setRho(params, algo, rhoOrS);
        end
        it0 = DefaultInitializeIterants(algo, ds.b_pad);
        it = optsolve(problem, algo, it0, ds.kernel, ds.b_pad, params);
        if isstruct(it), x = real(it.x); else, x = real(it); end
        x = min(max(x, 0), 1);
        xc = x(ds.crop(1):ds.crop(2), ds.crop(3):ds.crop(4));
        ps = psnr(xc, evalin('base','I'));
        ss = ssim(xc, evalin('base','I'));
    catch
    end
end

function p = setGamma(p, problem, g)
    if strcmp(problem,'l1'), p.gammal1 = g; else, p.gammal2 = g; end
end
function p = setT(p, algo, t)
    switch algo
        case 'douglasrachfordprimal',     p.tprimaldr     = t;
        case 'douglasrachfordprimaldual', p.tprimaldualdr = t;
        case 'admm',                      p.tadmm         = t;
        case 'chambollepock',             p.tchamb        = t;
    end
end
function p = setRho(p, algo, r)
    switch algo
        case 'douglasrachfordprimal',     p.rhoprimaldr     = r;
        case 'douglasrachfordprimaldual', p.rhoprimaldualdr = r;
        case 'admm',                      p.rhoadmm         = r;
    end
end