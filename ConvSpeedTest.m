% Clear and set random seed
clear; clc; close all;
rng(8);

outdir = 'results_compare';
if ~exist(outdir, 'dir'); mkdir(outdir); end

% Load and prep image
I = im2double(rgb2gray(imread('testimages/cameraman.jpg')));
I = imresize(I, [256, 256]);
[M, N] = size(I);

% Build the fixed test problem
kernel = fspecial('gaussian', [15, 15], 3); %15x15 with std of 3
kernel = kernel / sum(kernel(:));

apply_blur = @(img, ker) real(ifft2(fft2(img) .* ...
    eigValsForPeriodicConvOp(ker, size(img,1), size(img,2))));

b = apply_blur(I, kernel);
b = imnoise(b, 'salt & pepper', 0.05);
imwrite(I,               fullfile(outdir, 'truth.png'));
imwrite(min(max(b,0),1), fullfile(outdir, 'blurred.png'));

% ||A||^2 (for CP stability check)
eigK  = eigValsForPeriodicConvOp(kernel, M, N);
eigD1 = eigValsForPeriodicConvOp([-1,1]', M, N);
eigD2 = eigValsForPeriodicConvOp([-1,1],  M, N);
Anorm2 = max(max(abs(eigK).^2 + abs(eigD1).^2 + abs(eigD2).^2));
fprintf('||A||^2 = %.4f\n', Anorm2);

% Parameters per algorithm
common = struct();
common.maxiter = 2000; % Enough to see convergence
common.gammal1 = 0.01;
common.gammal2 = 0.005; % Not used since its an L1 problem

% CP
cp_prod  = 0.9999;
cp_ratio = 0.99;
common.schamb = sqrt(cp_prod / (Anorm2 * cp_ratio));
common.tchamb = cp_ratio * common.schamb;

% Primal DR
common.tprimaldr   = 0.3;
common.rhoprimaldr = 1;

% Primal-Dual DR
common.tprimaldualdr   = 10;
common.rhoprimaldualdr = 0.8;

% ADMM
common.tadmm   = 10;
common.rhoadmm = 0.3;

% Initializations
init_cp           = struct('x', b, 'y', zeros(M,N,3), 'z', b);
init_primaldr     = struct('z1', b, 'z2', zeros(M,N,3), ...
                           'x', b, 'y', zeros(M,N,3), ...
                           'b', zeros(M,N,3), 'u', b, 'v', zeros(M,N,3));
init_primaldualdr = struct('p', b, 'q', zeros(M,N,3), ...
                           'x', b, 'z', zeros(M,N,3), ...
                           'w', b, 'v', zeros(M,N,3));
init_admm         = struct('u', b, 'y', zeros(M,N,3), ...
                           'w', zeros(M,N), 'z', zeros(M,N,3), 'x', b);

% Algorithm registry
algos = {
    'chambollepock',              'Chambolle-Pock',         init_cp;
    'douglasrachfordprimal',      'Primal DR',              init_primaldr;
    'douglasrachfordprimaldual',  'Primal-Dual DR',         init_primaldualdr;
    'admm',                       'ADMM',                   init_admm;
};


%  Convergence comparison
%  - One run per algorithm with conv=true, measuring against I.
%  - Then one more run with conv=false to get the final image
%    for PSNR/SSIM reporting.

fprintf('\nConvergence comparison (%d iters each)\n', common.maxiter);

n_algos = size(algos,1);
results(n_algos) = struct('tag','','label','','dists',[], ...
                          'total_time',0,'time_per_iter',0, ...
                          'final_psnr',0,'final_ssim',0);

for a = 1:n_algos
    tag = algos{a,1}; label = algos{a,2}; init = algos{a,3};
    fprintf('\n--- %s ---\n', label);

    p = common;

    % Run 1: convergence trajectory against the original image I
    tic;
    dists = optsolve('l1', tag, init, kernel, b, p, I, true);
    time_conv = toc;

    % Run 2: get the final iterate for PSNR/SSIM
    tic;
    final = optsolve('l1', tag, init, kernel, b, p, I, false);
    time_image = toc;

    if isstruct(final), rec = real(final.x); else, rec = real(final); end
    rec = min(max(rec, 0), 1);

    total_time    = time_conv + time_image;
    time_per_iter = total_time / (2 * common.maxiter);

    psnr_val = psnr(rec, I);
    ssim_val = ssim(rec, I);

    results(a).tag           = tag;
    results(a).label         = label;
    results(a).dists         = dists;
    results(a).total_time    = total_time;
    results(a).time_per_iter = time_per_iter;
    results(a).final_psnr    = psnr_val;
    results(a).final_ssim    = ssim_val;

    fprintf('  wall-clock: %.2f s  (%.1f ms/iter)\n', total_time, time_per_iter*1000);
    fprintf('  final: PSNR %.2f dB,  SSIM %.3f\n', psnr_val, ssim_val);
    fprintf('  dist to I: start %.3e -> end %.3e\n', dists(1), dists(end));
end


%  Plotting
fig = figure('Name','Algorithm comparison','Position',[80 80 1400 500]);
cmap = lines(n_algos);
startIter = 5;  % skip the spike

% Floor subtraction: gives clean straight-line sublinear descent on log-log
d_floor  = min(arrayfun(@(r) min(r.dists), results));
eps_plot = 1e-6;   % avoid log(0)

% Panel (a): distance vs iteration  (log-log to reveal sublinear rate)
subplot(1,2,1); hold on;
for a = 1:n_algos
    k = startIter:length(results(a).dists);
    loglog(k, results(a).dists(startIter:end) - d_floor + eps_plot, ...
           'Color', cmap(a,:), 'LineWidth', 1.6, ...
           'DisplayName', results(a).label);
end
hold off;
set(gca,'XScale','log','YScale','log'); grid on;
xlabel('iteration k  (log scale)');
ylabel('||x^k - I||_2 - floor  (log scale)');
title('(a) Sublinear convergence vs iteration');
legend('Location','best');

% Panel (b): distance vs wall-clock time
subplot(1,2,2); hold on;
for a = 1:n_algos
    t_axis = (startIter:length(results(a).dists)) * results(a).time_per_iter;
    loglog(t_axis, results(a).dists(startIter:end) - d_floor + eps_plot, ...
           'Color', cmap(a,:), 'LineWidth', 1.6, ...
           'DisplayName', results(a).label);
end
hold off;
set(gca,'XScale','log','YScale','log'); grid on;
xlabel('wall-clock time (s)  (log scale)');
ylabel('||x^k - I||_2 - floor  (log scale)');
title('(b) Sublinear convergence vs wall-clock time');
legend('Location','best');

sgtitle(sprintf('Four-algorithm sublinear convergence comparison  (Gaussian blur + S&P 0.05, L1, %d iters)', common.maxiter));
saveas(fig, fullfile(outdir, 'algorithm_comparison.png'));

%  STEP 4 — Summary table
fprintf('\nSUMMARY\n');
fprintf('%-20s %-12s %-12s %-10s %-8s %-8s\n', ...
        'Algorithm','total time','ms/iter','final dist','PSNR','SSIM');
fprintf('%s\n', repmat('-',1,75));
for a = 1:n_algos
    fprintf('%-20s %-12.2f %-12.1f %-10.3e %-8.2f %-8.3f\n', ...
            results(a).label, ...
            results(a).total_time, ...
            results(a).time_per_iter*1000, ...
            results(a).dists(end), ...
            results(a).final_psnr, ...
            results(a).final_ssim);
end

fprintf('\nAll outputs in: %s/\n', outdir);