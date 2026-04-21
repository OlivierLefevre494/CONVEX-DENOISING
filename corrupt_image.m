% produces a 1x4 plot comparing the original image to three corrupted versions

clear; close all; clc;

% figure formatting
set(groot, 'defaultAxesFontSize', 18);
set(groot, 'defaultTextFontSize', 18);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultFigureColor', [1 1 1]);

imagePath = "testimages/cameraman.jpg";
I = PreprocessImage(imagePath);

%% Corruption settings
% Each entry specifies one blur/noise combination
corruptions(1).title = "Strong blur + mild noise";
corruptions(1).kernel = fspecial("gaussian", [9, 9], 8);
corruptions(1).noiseType = "salt & pepper";
corruptions(1).noiseArgs = {0.01};

corruptions(2).title = "Moderate blur + moderate noise";
corruptions(2).kernel = fspecial("gaussian", [9, 9], 4);
corruptions(2).noiseType = "salt & pepper";
corruptions(2).noiseArgs = {0.10};

corruptions(3).title = "Mild blur + strong noise";
corruptions(3).kernel = fspecial("gaussian", [9, 9], 0.5);
corruptions(3).noiseType = "salt & pepper";
corruptions(3).noiseArgs = {0.15};

% the first panel shows the original image, and the remaining panels show
% the corrupted images corresponding to the settings above.
fig = figure("Color", "w", "Position", [100, 100, 1800, 450]);
tiledlayout(1, numel(corruptions) + 1, "TileSpacing", "tight", "Padding", "compact");

nexttile;
imshow(I, []);
axis image off;
title("Original", "FontWeight", "bold", "Color", "k", "FontSize", 15);

for j = 1:numel(corruptions)
    % apply blur then noise
    blurred = imfilter(I, corruptions(j).kernel, "conv", "same");
    corrupted = imnoise(blurred, corruptions(j).noiseType, corruptions(j).noiseArgs{:});

    % clip intensities so the displayed image remains in [0,1].
    corrupted = min(max(corrupted, 0), 1);

    nexttile;
    imshow(corrupted, []);
    axis image off;
    title(corruptions(j).title, "FontWeight", "bold", "Color", "k", "FontSize", 15);
end
