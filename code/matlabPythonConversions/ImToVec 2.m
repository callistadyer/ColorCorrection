function [xVec] = ImToVec(xImage)
% Converts an image into one long vector using MATLAB ordering
%
% Syntax:
%    [xVec] = ImageToVec(xImage)
%
% Inputs:
%   xImage:   Input image of size H x W x C
%
% Outputs:
%   xVec:     Output vector of size (H*W*C) x 1
%
%   This ordering is the reference ordering that should also be matched in
%   Python when moving images between MATLAB and Python
%
% History:
%   03/10/2026  cmd    Wrote it
%
% Examples:
%{
xImage = rand(4,5,3);
xVec = ImageToVec(xImage);
%}

xVec = xImage(:);

end