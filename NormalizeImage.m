function [ NormalizeImage ] = NormalizeImage(image)
% Normalized the image sending pixel values to [0,1]

big = (image>1);
small = (image<0);
image(big) = 1.0;
image(small) = 0.0;
NormalizeImage = image; % Initialize the output variable
end
