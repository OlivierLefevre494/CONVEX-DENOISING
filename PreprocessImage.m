function [ preprocessimage ] = PreprocessImage(filename)
I = imread(filename); %% import image -> 2 gray -> double -> normalize 0-1
I = im2gray(I);
I = double(I(:, :, 1));
mn = min(I(:));
I = I-mn;
mx = max(I(:));
I = I/mx;
preprocessimage = I;

end