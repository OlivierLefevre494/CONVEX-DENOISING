function [ NormalizeImage ] = NormalizeImage(image)
% Normalized the image sending pixel values to [0,1]
[numRows, numCols] = size(image);

for i = 1:numRows
    for j = 1:numCols
        if image(i,j)>1
            image(i,j)=1.0;
        elseif image(i,j)<0
            image(i,j)=0.0;
        end
    end
end
NormalizeImage = image; % Initialize the output variable
end
