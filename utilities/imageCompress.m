function [Xout,w,h] = imageCompress(X,width,height, scale)
% X: the image collection, each column represents an image
% scale: the scale of image after downsampling
imgTemp = reshape(X(:,1),height,width);
imgTemp2 = imresize(imgTemp,scale);
[h,w] = size(imgTemp2);
Xout = zeros(h*w,size(X,2));
for i = 1:size(X,2)
    imgTemp = reshape(X(:,i),height,width);
    imgTemp2 = imresize(imgTemp,scale);
    Xout(:,i) = imgTemp2(:);
end

