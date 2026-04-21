%% USER INPUTS

clear

% path to input image
mypath = "testimages/cameraman.jpg";

% kernel and noise setting
kernel = fspecial('gaussian',9,4.0);
noisename = 'salt & pepper';
noisedens = 0.1;

% specify padding size for anti-aliasing
P=10; % THIS IS THE OPTIMAL PADDING FOR 9x9 KERNEL, CHANGE IT IF KERNEL SIZE CHANGES

% specify algorithm, problem and associated parameters
problem = 'l1';
algorithm = 'douglasrachfordprimaldual';
params = DefaultParams(); % can change this to this for user-defined params

%% READ, CORRUPT AND PADD IMAGE

myimage = PreprocessImage(mypath); % read image
blurredimage = imfilter(myimage, kernel); % blur image
blurredimage = imnoise(blurredimage, noisename, noisedens); % add noise

% we padd (replicate BCs) the corrupted image to mitigate aliasing errors
blurredimage = padarray(blurredimage, [P, P], "replicate");

%% RUN DENOISING ALGORITHM

iterants = DefaultInitializeIterants(algorithm, blurredimage); % can change the starting point if wanted...
outputimage = optsolve(problem, algorithm, iterants, kernel, blurredimage, params); % solve!

%% SHOW RESULTS

% chop padding
outputimage = outputimage(P+1:end-P, P+1:end-P);
blurredimage = blurredimage(P+1:end-P, P+1:end-P);

% set plotting parameters
set(groot,'defaultAxesFontSize',20);
set(groot,'defaultTextFontSize',20);
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelinterpreter','latex');
set(groot,'defaultLegendinterpreter','latex');
set(groot,'defaultLegendFontSize',20);
set(groot,'defaultfigurecolor',[1 1 1]);
set(0, 'DefaultLineLineWidth', 1);

figure(1)
tiledlayout(1, 3, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile; imshow(myimage,[]); title("Original Image")
nexttile; imshow(blurredimage,[]); title("Corrupted Image")
nexttile; imshow(outputimage,[]); title("Ouput")

