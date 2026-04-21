%% Set plotting params

set(groot,'defaultAxesFontSize',20);
set(groot,'defaultTextFontSize',20);
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelinterpreter','latex');
set(groot,'defaultLegendinterpreter','latex');
set(groot,'defaultLegendFontSize',20);
set(groot, 'DefaultTextColor', 'black');
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

kernel1 =fspecial('gaussian',30,3.01);
kernel2 = fspecial('motion', 18.4, 30);
noise1 = @(testimage) imnoise(testimage, 'salt & pepper',0.11);


%% List of possible boundary conditions
algorithm = ["chambollepock", "douglasrachfordprimal", "douglasrachfordprimaldual", "admm"];
boundary_conditions = ["None", "0", "Constant", "Replicate", "Symmetric", "Circular"];
padding_amount = [10, 20, 30, 40, 50, 60]; %% impact on speed?

function [ out ] = BoundaryCondition(boundary_condition, padding_amount, image, kernel, noise, algorithm)
CONST_VALUE = 0.5;

blurredimage = imfilter(image, kernel);
blurredimage = noise(blurredimage);
if strcmp(boundary_condition, "0")
    boundary_condition = 0;
end
if strcmp(boundary_condition, "None")
    padded_image=blurredimage;
elseif strcmp(boundary_condition, "Constant")
    padded_image = padarray(blurredimage, [padding_amount, padding_amount], CONST_VALUE, "both");
else
    padded_image = padarray(blurredimage, [padding_amount, padding_amount], lower(boundary_condition), "both");
end
iterants = DefaultInitializeIterants(algorithm, padded_image);
params = DefaultParams();
if strcmp(algorithm, "chambollepock")
    out = optsolve("l1", algorithm, iterants, kernel, padded_image, params, padded_image,false).x;
else
    out = optsolve("l1", algorithm, iterants, kernel, padded_image, params, padded_image,false);
end
% un-pad the image after processing
if strcmp(boundary_condition, "None")
    out = out;
else 
    out = out(padding_amount(1)+1:end-padding_amount(1), padding_amount(1)+1:end-padding_amount(1));
end
end


%% set up displayed figure
figure(1)
tiledlayout(2, 3, 'Padding', 'none', 'TileSpacing', 'compact');

% 
% %% run all boundary conditions
for i = 1:length(padding_amount)
     padamt = padding_amount(i);
     corruptblur1 = BoundaryCondition("Replicate", padamt, testimage, kernel1, noise1, "admm");
     MSEblur1 = immse(corruptblur1, testimage); % Calculate Mean Squared Error
     nexttile; imshow(corruptblur1,[]); title(sprintf(padamt+ "px, MSE = %.4f", MSEblur1))
end


%% for each combination of kernel/algorithm/padding type, run experiment, return table/struct of results
% results = table(); % Initialize an empty table to store results
% results2 = table();
% for i = 1:length(boundary_conditions)
%     for j = 1:length(padding_amount)
%         padding_amt = padding_amount(j);
%         for k = 1:length(algorithm)
%             for l = 1:2
%                 if l==1
%                     kernel = kernel1;
%                 else
%                     kernel = kernel2;
%                 end
%                 noise = noise1; 
%                 corruptImage = BoundaryCondition(boundary_conditions(i), padding_amt, testimage, kernel, noise, algorithm(k));
%                 MSE = immse(corruptImage, testimage); % Calculate Mean Squared Error
%                 if l==1
%                     results = [results; table(boundary_conditions(i), algorithm(k), MSE)]; % Store results
%                 else
%                     results2 = [results2; table(boundary_conditions(i), algorithm(k), MSE)]; % Store results
%                 end
%             end
%         end
%     end
% end
% 
% 
% 
% % Display results in a formatted table
% disp('Results for each combination of boundary conditions and algorithms:');
% disp(results);
% disp(results2);
