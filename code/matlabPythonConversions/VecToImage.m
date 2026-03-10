function [xImage] = VecToIm(xVec, imageSize)
% Converts a vector back into an image
%
% Syntax:
%    [xImage] = VecToImage(xVec, imageSize)
%
% Inputs:
%   xVec:        Input vector of size (H*W*C) x 1
%   imageSize:   Desired output image size [H W C]
%
% Outputs:
%   xImage:      Output image of size H x W x C
%
% Description:
%   This function reshapes a vector back into an image using MATLAB's
%   column-major ordering
%
%   It inverts ImageToVec, provided that imageSize matches the
%   size of the original image
%
%   That is, if
%       xVec = ImageToVec(xImage)
%   then
%       xImage = VecToImage(xVec, size(xImage))
%
% History:
%   03/10/2026  cmd    Wrote it.
%
% Examples:
%{
xImage = rand(4,5,3);
xVec = ImageToVec(xImage);
xImage2 = VecToImage(xVec, size(xImage));
isequal(xImage, xImage2)
%}

xImage = reshape(xVec, imageSize);

end