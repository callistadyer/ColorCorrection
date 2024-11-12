function [scaleFactor,lmsImageCalFormatOut,inGamut] = findDichromMappingScaleFactor(whichCone,lmsImageCalFormat,Disp,bScale)
% Take an LMS image in in cal format.  Find a scale factor to apply to the
% whichCone'th plane, so that after applying that scale factor, the LMS
% image is in gamut.  Return the scale factor, the scaled image, and
% whether the result is in fact in gamut.

% Initialize
scaleFactor = 1;
inGamut = false;
nTries = 100;
delta = scaleFactor/nTries;
lmsImageCalFormatOut = lmsImageCalFormat;

for i = 1:nTries
    % Scale with current scale factor
    lmsImageCalFormatOut(whichCone,:) = scaleFactor*lmsImageCalFormat(whichCone,:);

    % Is it in gamut.  Happiness if so.  If not, reduce scale factor and
    % loop.
    inGamut = checkGamut(lmsImageCalFormatOut,Disp,bScale);
    if (inGamut)
        break;
    else
        scaleFactor = scaleFactor-delta;
    end
end

% Let's barf here if it didn't work
if (~inGamut)
    error('Unable to find in gamut scale factor values');
end

end