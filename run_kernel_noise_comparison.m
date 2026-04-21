% Clear and set random seed
clear; clc; close all;
rng(8);

outdir = 'results_kernel_noise';
if ~exist(outdir, 'dir'); mkdir(outdir); end

% Load image
I = im2double(rgb2gray(imread('testimages/cameraman.jpg')));
I = imresize(I, [256, 256]);
[M, N] = size(I);
fprintf('Image: %dx%d\n', M, N);

% Kernels with blur helper
gauss_kernel  = fspecial('gaussian', [15, 15], 3);
gauss_kernel  = gauss_kernel / sum(gauss_kernel(:));
motion_kernel = fspecial('motion', 15, 30);
motion_kernel = motion_kernel / sum(motion_kernel(:));

apply_blur = @(img, ker) real(ifft2(fft2(img) .* ...
    eigValsForPeriodicConvOp(ker, size(img,1), size(img,2))));

% ||A||^2 per kernel
eigD1 = eigValsForPeriodicConvOp([-1,1]', M, N);
eigD2 = eigValsForPeriodicConvOp([-1,1],  M, N);

eigK_g   = eigValsForPeriodicConvOp(gauss_kernel, M, N);
Anorm2_g = max(max(abs(eigK_g).^2 + abs(eigD1).^2 + abs(eigD2).^2));

eigK_m   = eigValsForPeriodicConvOp(motion_kernel, M, N);
Anorm2_m = max(max(abs(eigK_m).^2 + abs(eigD1).^2 + abs(eigD2).^2));

% Init builders
make_init_cp   = @(b) struct('x', b, 'y', zeros(M,N,3), 'z', b);
make_init_pdr  = @(b) struct('z1', b, 'z2', zeros(M,N,3), ...
                             'x', b, 'y', zeros(M,N,3), ...
                             'b', zeros(M,N,3), 'u', b, 'v', zeros(M,N,3));
make_init_pddr = @(b) struct('p', b, 'q', zeros(M,N,3), ...
                             'x', b, 'z', zeros(M,N,3), ...
                             'w', b, 'v', zeros(M,N,3));
make_init_admm = @(b) struct('u', b, 'y', zeros(M,N,3), ...
                             'w', zeros(M,N), 'z', zeros(M,N,3), 'x', b);

% Test cases: kernel, noise, fidelity, gamma
cases = {
    'Gaussian blur + S&P noise',   gauss_kernel,  Anorm2_g, @(x) imnoise(x,'salt & pepper',0.05), 'l1', 0.01;
    'Gaussian blur + Gauss noise', gauss_kernel,  Anorm2_g, @(x) imnoise(x,'gaussian',0,0.01),    'l2', 0.009;
    'Motion blur + S&P noise',     motion_kernel, Anorm2_m, @(x) imnoise(x,'salt & pepper',0.05), 'l1', 0.01;
    'Motion blur + Gauss noise',   motion_kernel, Anorm2_m, @(x) imnoise(x,'gaussian',0,0.01),    'l2', 0.02;
};

algos      = {'chambollepock','douglasrachfordprimal','douglasrachfordprimaldual','admm'};
algo_labels = {'CP','Primal DR','PD-DR','ADMM'};
n_cases = size(cases, 1);
n_algos = length(algos);

psnr_table = zeros(n_cases, n_algos);
ssim_table = zeros(n_cases, n_algos);

% Run all combinations
for c = 1:n_cases
    label    = cases{c,1};
    kernel   = cases{c,2};
    Anorm2   = cases{c,3};
    noise_fn = cases{c,4};
    problem  = cases{c,5};
    gamma    = cases{c,6};

    fprintf('\nCase %d: %s [%s, gamma=%.4g]\n', c, label, problem, gamma);

    b = apply_blur(I, kernel);
    b = noise_fn(b);
    imwrite(min(max(b,0),1), fullfile(outdir, sprintf('case%d_blurred.png', c)));

    p = build_params(Anorm2, problem, gamma);

    for a = 1:n_algos
        algo = algos{a};

        switch algo
            case 'chambollepock',              init = make_init_cp(b);
            case 'douglasrachfordprimal',       init = make_init_pdr(b);
            case 'douglasrachfordprimaldual',   init = make_init_pddr(b);
            case 'admm',                        init = make_init_admm(b);
        end

        final = optsolve(problem, algo, init, kernel, b, p, I, false);
        if isstruct(final), rec = real(final.x); else, rec = real(final); end
        rec = min(max(rec, 0), 1);

        psnr_table(c, a) = psnr(rec, I);
        ssim_table(c, a) = ssim(rec, I);

        imwrite(rec, fullfile(outdir, sprintf('case%d_%s.png', c, algo)));
        fprintf('  %-12s  PSNR=%.2f  SSIM=%.3f\n', algo_labels{a}, psnr_table(c,a), ssim_table(c,a));
    end
end

% Summary tables
fprintf('\nPSNR TABLE (dB)\n');
fprintf('%-35s', 'Case');
for a = 1:n_algos, fprintf('%-12s', algo_labels{a}); end
fprintf('\n%s\n', repmat('-', 1, 35 + 12*n_algos));
for c = 1:n_cases
    fprintf('%-35s', cases{c,1});
    for a = 1:n_algos, fprintf('%-12.2f', psnr_table(c,a)); end
    fprintf('\n');
end

fprintf('\nSSIM TABLE\n');
fprintf('%-35s', 'Case');
for a = 1:n_algos, fprintf('%-12s', algo_labels{a}); end
fprintf('\n%s\n', repmat('-', 1, 35 + 12*n_algos));
for c = 1:n_cases
    fprintf('%-35s', cases{c,1});
    for a = 1:n_algos, fprintf('%-12.3f', ssim_table(c,a)); end
    fprintf('\n');
end

% Bar chart
figure('Name','Kernel/noise PSNR comparison','Position',[100 100 1000 500]);
bar(psnr_table);
set(gca, 'XTickLabel', {'Gauss+S&P (L1)','Gauss+GN (L2)','Motion+S&P (L1)','Motion+GN (L2)'});
ylabel('PSNR (dB)');
legend(algo_labels, 'Location', 'best');
title('PSNR across kernels and noise types (matched fidelity)');
grid on;
saveas(gcf, fullfile(outdir, 'kernel_noise_psnr_bar.png'));

% Before/after gallery
figure('Name','Before/after gallery','Position',[50 50 1400 900]);
for c = 1:n_cases
    [~, best_a] = max(psnr_table(c, :));

    subplot(n_cases, 3, (c-1)*3 + 1);
    imshow(imread(fullfile(outdir, sprintf('case%d_blurred.png', c))), []);
    if c == 1, title('Blurred + Noisy'); end
    ylabel(cases{c,1}, 'FontSize', 8, 'Interpreter', 'none');

    subplot(n_cases, 3, (c-1)*3 + 2);
    imshow(imread(fullfile(outdir, sprintf('case%d_%s.png', c, algos{best_a}))), []);
    if c == 1, title('Best Recovery'); end
    xlabel(sprintf('%s: %.2f dB', algo_labels{best_a}, psnr_table(c, best_a)));

    subplot(n_cases, 3, (c-1)*3 + 3);
    imshow(I, []);
    if c == 1, title('Ground Truth'); end
end
sgtitle('Kernel/noise comparison: blurred vs best recovery vs truth');
saveas(gcf, fullfile(outdir, 'before_after_gallery.png'));

fprintf('\nAll outputs in: %s/\n', outdir);

function p = build_params(Anorm2, problem, gamma)
    p = DefaultParams();
    p.maxiter = 2000;

    if strcmp(problem, 'l1')
        p.gammal1 = gamma;
    else
        p.gammal2 = gamma;
    end

    % CP step sizes
    cp_prod  = 0.9999;
    cp_ratio = 0.99;
    p.schamb = sqrt(cp_prod / (Anorm2 * cp_ratio));
    p.tchamb = cp_ratio * p.schamb;

    % DR variants
    p.tprimaldr     = 0.1;   p.rhoprimaldr     = 1.9;
    p.tprimaldualdr = 10;    p.rhoprimaldualdr = 1.9;

    % ADMM
    p.tadmm = 10.0;  p.rhoadmm = 1.35;
end