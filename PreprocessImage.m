function [ preprocessimage ] = PreprocessImage(filename)
I = imread(filename);
I = im2gray(I);
I = double(I(:, :, 1));
mn = min(I(:));
I = I-mn;
mx = max(I(:));
I = I/mx;
preprocessimage = I;

end