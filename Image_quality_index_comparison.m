%% Set plotting params

set(groot,'defaultAxesFontSize',20);
set(groot,'defaultTextFontSize',20);
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelinterpreter','latex');
set(groot,'defaultLegendinterpreter','latex');
set(groot,'defaultLegendFontSize',20);
set(groot,'defaultfigurecolor',[1 1 1]);
set(0, 'DefaultLineLineWidth', 1);


%% Why not MSE?

% Setup our testcase
testimage = PreprocessImage("testimages/cameraman.jpg");

% Corrupted image
kernel = fspecial('gaussian',9,4.0);
blurredimage = imfilter(testimage, kernel);
blurredimage = imnoise(blurredimage, 'salt & pepper', 0.1);
MSEref = mean((blurredimage - testimage).^2, 'all');
MSAref = mean(abs(blurredimage - testimage), 'all');
SSIMref = ssim(blurredimage, testimage);


%% Test with primal-dual results

algorithm = 'douglasrachfordprimaldual';
params.maxiter = 500;
params.tprimaldualdr = 3.0;
params.gammal1 = 0.25;
params.rhoprimaldualdr = 0.7;

iterants = DefaultInitializeIterants(algorithm, blurredimage);
problem = 'l1';
outputimage1 = optsolve(problem, algorithm, iterants, kernel, blurredimage, params);
MSE1 = mean((outputimage1 - testimage).^2, 'all')
MAE1 = mean(abs((outputimage1 - testimage)), 'all')
ssim1 = ssim(outputimage1, testimage)

algorithm = 'douglasrachfordprimaldual';
params.maxiter = 500;
params.tprimaldualdr = 0.8;
params.gammal1 = 0.049;
params.rhoprimaldualdr = 2.0;

iterants = DefaultInitializeIterants(algorithm, blurredimage);
problem = 'l1';
outputimage2 = optsolve(problem, algorithm, iterants, kernel, blurredimage, params);
MSE2 = mean((outputimage2 - testimage).^2, 'all')
MAE2 = mean(abs((outputimage2 - testimage)), 'all')
ssim2 = ssim(outputimage2, testimage)

figure(1)
tiledlayout(1, 3, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile; imshow(blurredimage,[]); title(sprintf("MSE = %.4f,    MAE = %.4f", MSEref, MSAref))
nexttile; imshow(outputimage1,[]); title(sprintf("MSE = %.4f,    MAE = %.4f", MSE1, MAE1))
nexttile; imshow(outputimage2,[]); title(sprintf("MSE = %.4f,    MAE = %.4f", MSE2, MAE2))

exportgraphics(gcf, "fig_PDDRS_imqual.png", "Resolution",500)
%% Test identical MSE
rng(1);

corruptblur1 = imfilter(testimage, fspecial('gaussian',50,4.5));
corruptsandp1 = imnoise(testimage, 'salt & pepper', 0.04);
corruptrot1 = imfilter(testimage, fspecial('motion',18.4,30));

MSEblur1 = mean((corruptblur1- testimage).^2, 'all');
MSEnoise1 = mean((corruptsandp1- testimage).^2, 'all');
MSErot1 = mean((corruptrot1 - testimage).^2, 'all');
ssimblur1 = ssim(corruptblur1, testimage);
ssimnoise1 = ssim(corruptsandp1, testimage);
ssimbrot1 = ssim(corruptrot1, testimage);

corruptblur2 = imfilter(testimage, fspecial('gaussian',50,3.01));
corruptsandp2 = imnoise(testimage, 'salt & pepper', 0.11);
corruptrot2 = imfilter(testimage, fspecial('motion',13.25,30));

MAEblur2= mean(abs(corruptblur2- testimage), 'all');
MAEnoise2 = mean(abs(corruptsandp2- testimage), 'all');
MAErot2 = mean(abs(corruptrot2 - testimage), 'all');
ssimblur2 = ssim(corruptblur2, testimage);
ssimnoise2 = ssim(corruptsandp2, testimage);
ssimbrot2 = ssim(corruptrot2, testimage);

figure(2)
tiledlayout(2, 3, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile; imshow(corruptblur1,[]); title(sprintf("MSE = %.4f", MSEblur1))
nexttile; imshow(corruptsandp1,[]); title(sprintf("MSE = %.4f", MSEnoise1))
nexttile; imshow(corruptrot1,[]); title(sprintf("MSE = %.4f", MSErot1))
nexttile; imshow(corruptblur2,[]); title(sprintf("MAE = %.4f", MAEblur2))
nexttile; imshow(corruptsandp2,[]); title(sprintf("MAE = %.4f", MAEnoise2))
nexttile; imshow(corruptrot2,[]); title(sprintf("MAE = %.4f", MAErot2))

%%
exportgraphics(gcf, "fig_corrupted_imqual.png", "Resolution",500)

%% FFT explanation
tiledlayout(1, 4, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile; imagesc(log(fftshift(abs(fft2(testimage)))), [-5,10]); axis equal; xlim([0 256]); ylim([0 256]); xticks([]); xlabel("$k_x$"); yticks([]); ylabel("$k_y$"); title("Reference"); colorbar
nexttile; imagesc(log(fftshift(abs(fft2(corruptblur1)))), [-5,10]); axis equal; xlim([0 256]); ylim([0 256]); xticks([]); xlabel("$k_x$"); yticks([]); ylabel("$k_y$"); title("Gaussian Blur"); colorbar
nexttile; imagesc(log(fftshift(abs(fft2(corruptsandp1)))), [-5,10]); axis equal; xlim([0 256]); ylim([0 256]); xticks([]); xlabel("$k_x$"); yticks([]); ylabel("$k_y$"); title("S\&P Noise"); colorbar
nexttile; imagesc(log(fftshift(abs(fft2(corruptrot1)))), [-5,10]); axis equal; xlim([0 256]); ylim([0 256]); xticks([]); xlabel("$k_x$"); yticks([]); ylabel("$k_y$"); title("Motion Blur"); colorbar