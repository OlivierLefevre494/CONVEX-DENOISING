% Hyperparameter visualization for Section 4.1
% Produces three figures. Each figure varies one parameter while the
% other two are held fixed.

clear; close all; clc;

rng(1);

% figure formatting 
set(groot, 'defaultAxesFontSize', 20);
set(groot, 'defaultTextFontSize', 20);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultLegendFontSize', 20);
set(groot, 'defaultFigureColor', [1 1 1]);
set(groot, 'defaultLineLineWidth', 1);

algorithm = "admm";
problem = "l1";               
maxiter = 500;

% fixed parameter values 
fixedGamma = 0.03;
fixedRho = 0.5;
fixedT = 1.0;

% parameter values to vary 
gammaSweep = [0.01, 0.02, 0.03, 0.04];
rhoSweep = [1.50, 1.52, 1.54, 1.56];
tSweep = [0.5, 1.0, 1.5, 2.0];

imagePath = "testimages/cameraman.jpg";
kernel = fspecial("gaussian", [9, 9], 4);
noiseType = "salt & pepper";
noiseDensity = 0.2;

% option to use padding
usePadding = true;
pad = 10;

% corrupt image
I = PreprocessImage(imagePath);
b = imfilter(I, kernel, "conv", "same");
b = imnoise(b, noiseType, noiseDensity);

if usePadding
    B = padarray(b, [pad, pad], "replicate", "both");
else
    B = b;
end

% define the three parameter sweeps
sweeps(1).name = "$\gamma$";
sweeps(1).values = gammaSweep;
sweeps(1).gammaVals = gammaSweep;
sweeps(1).rhoVals = fixedRho * ones(size(gammaSweep));
sweeps(1).tVals = fixedT * ones(size(gammaSweep));
sweeps(1).fixedText = sprintf("$\\rho=%.2g$, $t=%.2g$", fixedRho, fixedT);

sweeps(2).name = "$\rho$";
sweeps(2).values = rhoSweep;
sweeps(2).gammaVals = fixedGamma * ones(size(rhoSweep));
sweeps(2).rhoVals = rhoSweep;
sweeps(2).tVals = fixedT * ones(size(rhoSweep));
sweeps(2).fixedText = sprintf("$\\gamma=%.3g$, $t=%.2g$", fixedGamma, fixedT);

sweeps(3).name = "$t$";
sweeps(3).values = tSweep;
sweeps(3).gammaVals = fixedGamma * ones(size(tSweep));
sweeps(3).rhoVals = fixedRho * ones(size(tSweep));
sweeps(3).tVals = tSweep;
sweeps(3).fixedText = sprintf("$\\gamma=%.3g$, $\\rho=%.2g$", fixedGamma, fixedRho);

% run sweeps
params = DefaultParams();
params.maxiter = maxiter;
numSweeps = numel(sweeps);
numValues = numel(gammaSweep);
reconstructions = cell(numSweeps, numValues);

for row = 1:numSweeps
    for col = 1:numValues
        gamma = sweeps(row).gammaVals(col);
        rho = sweeps(row).rhoVals(col);
        t = sweeps(row).tVals(col);

        params = setAlgorithmParams(params, algorithm, t, rho);
        if problem == "l1"
            params.gammal1 = gamma;
        else
            params.gammal2 = gamma;
        end

        fprintf("Running %s: %s, gamma = %.4g, rho = %.4g, t = %.4g\n", ...
            algorithm, problem, gamma, rho, t);

        iterants = DefaultInitializeIterants(algorithm, B);
        result = optsolve(problem, algorithm, iterants, kernel, B, params);
        result = real(result);

        if usePadding
            result = result(1+pad:end-pad, 1+pad:end-pad);
        end
        result = min(max(result, 0), 1);

        reconstructions{row, col} = result;
    end
end

% plots three figures, one for each sweep
for row = 1:numSweeps
    fig = figure("Color", "w", "Position", [100, 100, 1500, 450]);
    layout = tiledlayout(1, numValues, "TileSpacing", "tight", "Padding", "compact");

    for col = 1:numValues
        nexttile;
        imshow(reconstructions{row, col}, []);
        axis image off;

        title(sprintf("%s $= %.3g$", sweeps(row).name, sweeps(row).values(col)), ...
            "FontWeight", "bold", "Color", "k", "FontSize", 16);
    end
end


function params = setAlgorithmParams(params, algorithm, t, rho)
    switch algorithm
        case "douglasrachfordprimal"
            params.tprimaldr = t;
            params.rhoprimaldr = rho;
        case "douglasrachfordprimaldual"
            params.tprimaldualdr = t;
            params.rhoprimaldualdr = rho;
        case "admm"
            params.tadmm = t;
            params.rhoadmm = rho;
        case "chambollepock"
            params.tchamb = t;
            % Chambolle-Pock has no rho parameter in this implementation.
        otherwise
            error('Unknown algorithm "%s".', algorithm);
    end
end
