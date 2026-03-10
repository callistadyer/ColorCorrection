function [calFormat] = VecToCalFormat(xVec, calFormatSize)
% Converts a vector back into a CalFormat matrix 
%
% Syntax:
%    [calFormat] = VecToCalFormat(xVec, calFormatSize)
%
% Inputs:
%   xVec:            Input vector of size (3*N) x 1
%   calFormatSize:   Desired output CalFormat size [3 N]
%
% Outputs:
%   calFormat:       Output CalFormat matrix of size 3 x N
%
% Description:
%   This function reshapes a vector back into a CalFormat matrix using
%   MATLAB's native column-major ordering.
%
%   It inverts CalFormatToVec, provided that calFormatSize matches
%   the size of the original CalFormat matrix.
%
%   That is, if
%       xVec = CalFormatToVec(calFormat)
%   then
%       calFormat = VecToCalFormat(xVec, size(calFormat))
%
% History:
%   03/10/2026  cmd    Wrote it.
%
% Examples:
%{
calFormat = [1 2 3; 10 20 30; 100 200 300];
xVec = CalFormatToVec(calFormat);
calFormat2 = VecToCalFormat(xVec, size(calFormat));
isequal(calFormat, calFormat2)
%}

calFormat = reshape(xVec, calFormatSize);

end