%% Set plotting params

set(groot,'defaultAxesFontSize',20);
set(groot,'defaultTextFontSize',20);
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelinterpreter','latex');
set(groot,'defaultLegendinterpreter','latex');
set(groot,'defaultLegendFontSize',20);
set(groot,'defaultfigurecolor',[1 1 1]);
set(0, 'DefaultLineLineWidth', 1);

rng(1);



testimage = PreprocessImage("testimages/cameraman.jpg");


% corruptblur2 = imfilter(testimage, fspecial('gaussian',50,3.01));
% corruptsandp2 = imnoise(testimage, 'salt & pepper', 0.11);
% corruptrot2 = imfilter(testimage, fspecial('motion',13.25,30));
% corruptblur1 = imfilter(testimage, fspecial('gaussian',50,4.5));
% corruptsandp1 = imnoise(testimage, 'salt & pepper', 0.04);
% corruptrot1 = imfilter(testimage, fspecial('motion',18.4,30));

kernel =fspecial('gaussian',50,4.5);
noise = @(testimage) imnoise(testimage, 'salt & pepper',0.11);

%% List of possible boundary conditions

boundary_conditions = ["none", "0", "constant", "replicate", "symmetric", "circular"];
padding_amount = [20]; %% impact on speed?

function [ out ] = BoundaryCondition(boundary_condition, padding_amount, image, kernel, noise)
CONST_VALUE = 0.5;
if strcmp(boundary_condition, "0")
    boundary_condition = 0;
end
if strcmp(boundary_condition, "none")
    padded_image=image;
elseif strcmp(boundary_condition, "constant")
    padded_image = padarray(image, [padding_amount, padding_amount], CONST_VALUE, "both");
else
    padded_image = padarray(image, [padding_amount, padding_amount], boundary_condition, "both");
end
blurredimage = imfilter(padded_image, kernel);
blurredimage = noise(blurredimage);
iterants = DefaultInitializeIterants("douglasrachfordprimaldual", blurredimage);
params.maxiter = 500;
params.tprimaldualdr = 0.8;
params.gammal1 = 0.01;
params.rhoprimaldualdr = 1.0;
out = optsolve("l1", "douglasrachfordprimaldual", iterants, kernel, blurredimage, params, padded_image,false);
% un-pad the image after processing
if strcmp(boundary_condition, "none")
    out = out;
else 
    out = out(padding_amount(1)+1:end-padding_amount(1), padding_amount(1)+1:end-padding_amount(1));
end
end


%% set up displayed figure
figure(1)
tiledlayout(2, 3, 'Padding', 'none', 'TileSpacing', 'compact');


%% run all boundary conditions
for i = 1:length(boundary_conditions)
    current_condition = boundary_conditions(i);
    corruptblur1 = BoundaryCondition(current_condition, padding_amount, testimage, kernel, noise);
    MSEblur1 = immse(corruptblur1, testimage); % Calculate Mean Squared Error
    nexttile; imshow(corruptblur1,[]); title(sprintf(current_condition+ ", MSE = %.4f", MSEblur1))

end