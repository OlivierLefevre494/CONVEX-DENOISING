clear; clc; close all;
rng(42);

outdir = 'results_l1_vs_l2';
if ~exist(outdir, 'dir'); mkdir(outdir); end

% Load image
I = im2double(rgb2gray(imread('testimages/cameraman.jpg')));
I = imresize(I, [256, 256]);
[M, N] = size(I);
fprintf('Image: %dx%d\n', M, N);

% Kernel
kernel = fspecial('gaussian', [15, 15], 3);
kernel = kernel / sum(kernel(:));

apply_blur = @(img, ker) real(ifft2(fft2(img) .* ...
    eigValsForPeriodicConvOp(ker, size(img,1), size(img,2))));

% ||A||^2
eigK  = eigValsForPeriodicConvOp(kernel, M, N);
eigD1 = eigValsForPeriodicConvOp([-1,1]', M, N);
eigD2 = eigValsForPeriodicConvOp([-1,1],  M, N);
Anorm2 = max(max(abs(eigK).^2 + abs(eigD1).^2 + abs(eigD2).^2));

% Two noise scenarios
b_sp = apply_blur(I, kernel);
b_sp = imnoise(b_sp, 'salt & pepper', 0.05);

b_gn = apply_blur(I, kernel);
b_gn = imnoise(b_gn, 'gaussian', 0, 0.01);

imwrite(min(max(b_sp,0),1), fullfile(outdir, 'blurred_sp.png'));
imwrite(min(max(b_gn,0),1), fullfile(outdir, 'blurred_gn.png'));
imwrite(I,                  fullfile(outdir, 'truth.png'));

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

% Algorithms
algos = {'chambollepock', 'douglasrachfordprimal', 'douglasrachfordprimaldual', 'admm'};
algo_labels = {'CP', 'Primal DR', 'PD-DR', 'ADMM'};
n_algos = length(algos);

% Test grid
noise_cases = {
    'S&P 0.05',      b_sp;
    'Gaussian 0.01',  b_gn;
};
problems = {'l1', 'l2'};

n_noise = size(noise_cases, 1);
n_prob  = length(problems);

% Storage: psnr_results(noise, problem, algo)
psnr_results = zeros(n_noise, n_prob, n_algos);
ssim_results = zeros(n_noise, n_prob, n_algos);
images       = cell(n_noise, n_prob, n_algos);

for ni = 1:n_noise
    noise_label = noise_cases{ni, 1};
    b = noise_cases{ni, 2};

    for pi = 1:n_prob
        prob = problems{pi};
        fprintf('\n=== %s, %s ===\n', noise_label, upper(prob));

        p = build_params(Anorm2, b);

        for a = 1:n_algos
            algo = algos{a};

            switch algo
                case 'chambollepock',              init = make_init_cp(b);
                case 'douglasrachfordprimal',       init = make_init_pdr(b);
                case 'douglasrachfordprimaldual',   init = make_init_pddr(b);
                case 'admm',                        init = make_init_admm(b);
            end

            final = optsolve(prob, algo, init, kernel, b, p, I, false);
            if isstruct(final), rec = real(final.x); else, rec = real(final); end
            rec = min(max(rec, 0), 1);

            psnr_results(ni, pi, a) = psnr(rec, I);
            ssim_results(ni, pi, a) = ssim(rec, I);
            images{ni, pi, a} = rec;

            fname = sprintf('%s_%s_%s.png', noise_label, prob, algo);
            fname = regexprep(fname, '[^a-zA-Z0-9_.]', '_');
            imwrite(rec, fullfile(outdir, fname));

            fprintf('  %-12s  PSNR=%.2f  SSIM=%.3f\n', ...
                    algo_labels{a}, psnr_results(ni,pi,a), ssim_results(ni,pi,a));
        end
    end
end

%Summary tables
fprintf('\nPSNR TABLE (dB)\n');
for ni = 1:n_noise
    fprintf('\n--- %s ---\n', noise_cases{ni,1});
    fprintf('%-12s', 'Algorithm');
    for pi = 1:n_prob, fprintf('%-12s', upper(problems{pi})); end
    fprintf('\n%s\n', repmat('-', 1, 12 + 12*n_prob));
    for a = 1:n_algos
        fprintf('%-12s', algo_labels{a});
        for pi = 1:n_prob, fprintf('%-12.2f', psnr_results(ni,pi,a)); end
        fprintf('\n');
    end
end

fprintf('\nSSIM TABLE\n');
for ni = 1:n_noise
    fprintf('\n--- %s ---\n', noise_cases{ni,1});
    fprintf('%-12s', 'Algorithm');
    for pi = 1:n_prob, fprintf('%-12s', upper(problems{pi})); end
    fprintf('\n%s\n', repmat('-', 1, 12 + 12*n_prob));
    for a = 1:n_algos
        fprintf('%-12s', algo_labels{a});
        for pi = 1:n_prob, fprintf('%-12.3f', ssim_results(ni,pi,a)); end
        fprintf('\n');
    end
end

% Bar chart: grouped by noise, L1 vs L2 side by side
% For each noise type, pick the best algorithm's PSNR for L1 and L2
figure('Name','L1 vs L2 PSNR','Position',[100 100 900 500]);

% Reshape: for each noise x problem, take the max across algorithms
best_psnr = zeros(n_noise, n_prob);
best_algo = cell(n_noise, n_prob);
for ni = 1:n_noise
    for pi = 1:n_prob
        [best_psnr(ni,pi), idx] = max(psnr_results(ni,pi,:));
        best_algo{ni,pi} = algo_labels{idx};
    end
end

bar(best_psnr);
set(gca, 'XTickLabel', {noise_cases{1,1}, noise_cases{2,1}});
ylabel('PSNR (dB)');
legend({'L1', 'L2'}, 'Location', 'best');
title('Best PSNR: L1 vs L2 formulation by noise type');
grid on;

% Annotate which algorithm won
for ni = 1:n_noise
    for pi = 1:n_prob
        x_pos = ni + (pi - 1.5) * 0.28;
        text(x_pos, best_psnr(ni,pi) + 0.3, best_algo{ni,pi}, ...
             'HorizontalAlignment', 'center', 'FontSize', 8);
    end
end

saveas(gcf, fullfile(outdir, 'l1_vs_l2_best_psnr.png'));

% Full bar chart: all algorithms x L1/L2 x noise
figure('Name','L1 vs L2 all algorithms','Position',[50 50 1400 500]);

for ni = 1:n_noise
    subplot(1, 2, ni);
    data = squeeze(psnr_results(ni, :, :))';  % n_algos x n_prob
    bar(data);
    set(gca, 'XTickLabel', algo_labels);
    ylabel('PSNR (dB)');
    legend({'L1', 'L2'}, 'Location', 'best');
    title(sprintf('%s', noise_cases{ni,1}));
    grid on;
end
sgtitle('L1 vs L2: PSNR by algorithm and noise type');
saveas(gcf, fullfile(outdir, 'l1_vs_l2_all_algos.png'));

% 2x2 image comparison (best L1 vs best L2 per noise)
figure('Name','L1 vs L2 images','Position',[100 100 900 900]);
for ni = 1:n_noise
    for pi = 1:n_prob
        [~, best_a] = max(psnr_results(ni, pi, :));
        subplot(2, 2, (ni-1)*2 + pi);
        imshow(images{ni, pi, best_a}, []);
        title(sprintf('%s, %s\n%s: PSNR %.2f', ...
              noise_cases{ni,1}, upper(problems{pi}), ...
              best_algo{ni,pi}, best_psnr(ni,pi)));
    end
end
sgtitle('Best recovery: L1 vs L2 across noise types');
saveas(gcf, fullfile(outdir, 'l1_vs_l2_images.png'));

fprintf('\nAll outputs in: %s/\n', outdir);

function p = build_params(Anorm2, b)
    p = DefaultParams();
    p.maxiter = 1000;
    p.gammal1 = 0.01;
    p.gammal2 = 0.005;
    
    % CP
    cp_prod  = 0.9999;
    cp_ratio = 0.99;
    p.schamb = sqrt(cp_prod / (Anorm2 * cp_ratio));
    p.tchamb = cp_ratio * p.schamb;
    
    % Primal DR
    p.tprimaldr   = 0.1;
    p.rhoprimaldr = 1.9;
    
    % Primal-Dual DR
    p.tprimaldualdr   = 10;
    p.rhoprimaldualdr = 1.9;
    
    % ADMM
    p.tadmm   = 10.0;
    p.rhoadmm = 1.35;
end