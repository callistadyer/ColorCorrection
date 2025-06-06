function [scaleFactor,lmsImageCalFormatOut,inGamut] = findDichromMappingScaleFactor(whichCone,lmsImageCalFormat,Disp,bScale)
% Take an LMS image in in cal format.  Find a scale factor to apply to the
% whichCone'th plane, so that after applying that scale factor, the LMS
% image is in gamut.  Return the scale factor, the scaled image, and
% whether the result is in fact in gamut.

% Initialize
scaleFactor = 1;
inGamut = false;
nTries = 300;
delta = scaleFactor / nTries;
lmsImageCalFormatOut = lmsImageCalFormat;

% Track the maximum value of rgbImageCalFormat
maxValues = zeros(1, nTries);

for i = 1:nTries
    % Scale with current scale factor
    lmsImageCalFormatOut(whichCone, :) = scaleFactor * lmsImageCalFormat(whichCone, :);

    % Check if it is in gamut
    [inGamut, rgbImageCalFormat] = checkGamut(lmsImageCalFormatOut, Disp, bScale);

    % Store the maximum value of rgbImageCalFormat
    maxValues(i) = max(rgbImageCalFormat(:));

    % If in gamut, exit
    if inGamut
        break;
    else
        scaleFactor = scaleFactor - delta;
    end

    % Check if max value at 100th iteration is larger than the first
    if i == 100
        if maxValues(100) > maxValues(1)
            % Restart the scaling process with reversed direction
            scaleFactor = 1; % Reset scale factor
            delta = -delta; % Reverse scaling direction
            maxValues = zeros(1, nTries); % Reset maxValues tracking
            i = 0; % Restart the loop
        end
    end
end

% Check if it succeeded
if ~inGamut
    error('Unable to find in-gamut scale factor values');
end

end


% % Initialize
% scaleFactor = 1;
% inGamut = false;
% nTries = 300;
% delta = scaleFactor/nTries;
% lmsImageCalFormatOut = lmsImageCalFormat;
% 
% for i = 1:nTries
%     % Scale with current scale factor
%     lmsImageCalFormatOut(whichCone,:) = scaleFactor*lmsImageCalFormat(whichCone,:);
% 
%     % Is it in gamut.  Happiness if so.  If not, reduce scale factor and
%     % loop.
%     [inGamut,rgbImageCalFormat] = checkGamut(lmsImageCalFormatOut,Disp,bScale);
%     if (inGamut)
%         break;
%     else
%         scaleFactor = scaleFactor-delta;
%     end
% end
% 
% % Let's barf here if it didn't work
% if (~inGamut)
%     error('Unable to find in gamut scale factor values');
% end
% 
% end