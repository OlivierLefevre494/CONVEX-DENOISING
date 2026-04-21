% Chambolle--Pock step-size heatmap.
% This script sweeps over a grid of dual/primal step sizes (s,t), runs the
% Chambolle--Pock solver for each pair, and records the final L2
% reconstruction error. The black curve is the theoretical stability
% boundary s*t*||A||^2 = 1.

clear; close all; clc;

% figure formatting
set(groot, 'defaultAxesFontSize', 18);
set(groot, 'defaultTextFontSize', 18);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultFigureColor', [1 1 1]);

% choose the solver and reconstruction model.
algorithm = "chambollepock";
problem = "l1";
gamma = 0.03;

% number of Chambolle--Pock iterations performed at each grid point
% we found that >500 produced a more meaningful visualization
maxiter = 1000;

% grid of step sizes to be plotted in (s,t)-space
sGrid = linspace(0.03, 0.50, 15);
tGrid = linspace(0.10, 2.00, 15);

% test image and corruption model 
imagePath = "testimages/cameraman.jpg";
kernel = fspecial("gaussian", [9, 9], 4);
noiseType = "salt & pepper";
noiseDensity = 0.10;

% option to include padding
usePadding = true;
pad = 9;

% load the clean image, blur it, and add noise to form the observed data b.
rng(1);
I = PreprocessImage(imagePath);
b = imfilter(I, kernel, "conv", "same");
b = imnoise(b, noiseType, noiseDensity);

if usePadding
    B = padarray(b, [pad, pad], "replicate", "both");
else
    B = b;
end

% compute ||A||^2 for the stability boundary
[numRows, numCols] = size(B);
eigArry_K = eigValsForPeriodicConvOp(kernel, numRows, numCols);
eigArry_D1 = eigValsForPeriodicConvOp([-1, 1]', numRows, numCols);
eigArry_D2 = eigValsForPeriodicConvOp([-1, 1], numRows, numCols);
L2 = max(abs(eigArry_K).^2 + abs(eigArry_D1).^2 + abs(eigArry_D2).^2, [], "all");


params = DefaultParams();
params.maxiter = maxiter;
params.gammal1 = gamma;

% `errors(row,col)` stores the final reconstruction error at the step-size
% pair (sGrid(row), tGrid(col)). `stabilityCheck(row,col)` stores the
% quantity s*t*||A||^2 associated with the convergence condition.
errors = NaN(numel(sGrid), numel(tGrid));
stabilityCheck = NaN(numel(sGrid), numel(tGrid));

numTrials = numel(sGrid) * numel(tGrid);
trial = 1;

for row = 1:numel(sGrid)
    for col = 1:numel(tGrid)
        s = sGrid(row);
        t = tGrid(col);

        % update the Chambolle--Pock step sizes for this trial.
        params.schamb = s;
        params.tchamb = t;
        stabilityCheck(row, col) = s * t * L2;

        fprintf("Trial %3d/%3d: s = %.4g, t = %.4g, s*t*||A||^2 = %.4g\n", ...
            trial, numTrials, s, t, stabilityCheck(row, col));

        try
            % reinitialize the iterates for each run so all grid points are
            % compared from the same starting state.
            iterants = DefaultInitializeIterants("chambollepock", B);
            result = optsolve(problem, algorithm, iterants, kernel, B, params);

            % extract the reconstruction and undo the padding before
            % comparing against the original image.
            reconstruction = real(result.x);
            if usePadding
                reconstruction = reconstruction(1+pad:end-pad, 1+pad:end-pad);
            end
            reconstruction = min(max(reconstruction, 0), 1);

            % only record the error if the solver produced a finite image.
            if all(isfinite(reconstruction(:)))
                errors(row, col) = norm(reconstruction(:) - I(:), 2);
            end
        catch ME
            % keep the sweep running even if a particular step-size pair
            % fails or diverges numerically.
            warning("Failed at s = %.4g, t = %.4g: %s", s, t, ME.message);
        end

        trial = trial + 1;
    end
end

% plot the heatmap and the stability boundary
fig = figure("Color", "w", "Position", [100, 100, 900, 650]);
imagesc(tGrid, sGrid, errors);
set(gca, "YDir", "normal", "Color", "w", "XColor", "k", "YColor", "k");
colormap(parula);
cb = colorbar;
cb.Label.String = "L2 error";
cb.Label.Interpreter = "latex";
cb.Color = "k";
hold on;

sBoundary = linspace(min(sGrid), max(sGrid), 300);
tBoundary = 1 ./ (sBoundary .* L2);
valid = tBoundary >= min(tGrid) & tBoundary <= max(tGrid);
plot(tBoundary(valid), sBoundary(valid), "k-", "LineWidth", 2.5);

xlabel("$t$");
ylabel("$s$");
title("Chambolle--Pock Step Sizes", "FontWeight", "bold", "Color", "k");
box on;
