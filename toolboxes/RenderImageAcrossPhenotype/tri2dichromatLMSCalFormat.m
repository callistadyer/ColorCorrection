function [diLMSCalFormat] = tri2dichromatLMSCalFormat(triLMSCalFormat,renderType,Disp,bScale)

% function takes in trichromat lms values and converts them into dichromat lms values
%
% Syntax:
%   [diLMSCalFormat] = tri2dichromatLMSCalFormat(triLMSCalFormat,renderType,Disp,bScale)
%
% Description:
%
% Inputs:
%   triLMSCalFormat    LMS values of the image for a trichromat
%   renderType         Type of dichromacy
%                        Deuteranopia
%                        Protanopia
%                        Tritanopia
%
% Outputs:
%   diLMSCalFormat     LMS values for dichromat type specified by renderType
%
% Optional key/value pairs:
%   None

switch (renderType)
    case 'Deuteranopia' % m cone deficiency
        cbType = 2;
    case 'Protanopia'   % l cone deficiency
        cbType = 1;
    case 'Tritanopia'   % s cone deficiency
        cbType = 3;
end

% LMS --> Linear RGB (so we can go from RGB --> XYZ)
triRGBlinCalFormat = LMS2rgbLinCalFormat(triLMSCalFormat,Disp,0);

% Matrix to convert from rgb to xyz
% Note: this matrix must be applied on the LEFT!!
M_rgb2xyz = Disp.T_xyz*Disp.P_monitor;
M_xyz2rgb = inv(M_rgb2xyz);

% Linear RGB --> XYZ (Brettel takes in XYZ)
triXYZCalFormat = M_rgb2xyz * triRGBlinCalFormat;
% Cal Format --> Image Format
triXYZImgFormat = CalFormatToImage(triXYZCalFormat,Disp.m,Disp.n);

%%%%%%%%% Dichromat Simulation (Brettel) %%%%%%%%%%%%%%%%%%%%%%
[diXYZ] = DichromatSimulateBrettel(triXYZImgFormat, cbType, []);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Image Format --> CalFormat
diXYZCalFormat = ImageToCalFormat(diXYZ);
% XYZ --> Linear RGB
diRGBLinCalFormat = M_xyz2rgb * diXYZCalFormat;

% Quick snipping if the vals are only over by a small amount
cutoffUnder = 1.1;
if (max(diRGBLinCalFormat(:)) > 1) && (max(diRGBLinCalFormat(:)) < cutoffUnder)
    diRGBLinCalFormat(diRGBLinCalFormat>1) = .99; % For some reason, still goes out of gamut if I make these vals == 1, so .99 seems to work
end
% Linear RGB --> LMS 
diLMSCalFormat = rgbLin2LMSCalFormat(diRGBLinCalFormat,Disp,1,0);

% CHECK IF MODULATED LMS IS IN GAMUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inGamut = checkGamut(diLMSCalFormat,Disp,bScale);
if inGamut == 0
    error(['tri2dichromatLMSCalFormat: WARNING! rgb values are out of gamut... lmsImage_mod values outside of the range [0 1]']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














% This function doesn't work... my attempt at Brettel simulation (use DichromatSimulateBrettel instead)  
% diLMSCalFormat = simulateDichromatBrettel(lmsImage,renderType,Disp);

%%%% OLD WAY OF SIMULATING DICHROMAT... THIS PUSHES STUFF OUT OF GAMUT, AND IS ALSO KINDA ARBITRARY %%%% 
% NOTE WELL: (outdated)
% To get the renderings to look right, we are using the mean value of
% at the missing cone plane to get the scale of our substitution
% right.  Presumably this could be done once by analyzing a full ensemble
% of images, but we are going to worry about that later.
% l_cone = triLMSCalFormat(1,:);
% m_cone = triLMSCalFormat(2,:);
% s_cone = triLMSCalFormat(3,:);
% % Make dichromat manipulation - missing cone
% switch (renderType)
%     case 'Deuteranopia' % m cone deficiency
%         lmsImageCalFormat(2,:)       = l_cone; 
%         [scaleFactor,dichromatLMSCalFormat,inGamut] = findDichromMappingScaleFactor(2,lmsImageCalFormat,Disp,bScale);
% 
%     case 'Protanopia'   % l cone deficiency
%         % CHECK: this had a 2 in the index... i changed it to 1. correct?
%         lmsImageCalFormat(1,:)       = m_cone; 
%         [scaleFactor,dichromatLMSCalFormat,inGamut] = findDichromMappingScaleFactor(1,lmsImageCalFormat,Disp,bScale); 
% 
%     case 'Tritanopia'   % s cone deficiency
%         lmsImageCalFormat(3,:)       = (m_cone + l_cone)/2;
%         [scaleFactor,dichromatLMSCalFormat,inGamut] = findDichromMappingScaleFactor(3,lmsImageCalFormat,Disp,bScale);
% 
% end

%dichromatLMSCalFormat = lmsImageCalFormat;

end