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

l_cone = triLMSCalFormat(1,:);
m_cone = triLMSCalFormat(2,:);
s_cone = triLMSCalFormat(3,:);

% NOTE WELL:
% To get the renderings to look right, we are using the mean value of
% at the missing cone plane to get the scale of our substitution
% right.  Presumably this could be done once by analyzing a full ensemble
% of images, but we are going to worry about that later.

%%%%%%%%%%%%%%%%%%% CALLISTA WORK ON THIS!!!!
% disp(['callista work on this! choose how to simulate dichromat']);
lmsImage = CalFormatToImage(triLMSCalFormat,Disp.m,Disp.n);
triRGBlinCalFormat = LMS2rgbLinCalFormat(triLMSCalFormat,Disp,0);
triRGBlinImgFormat = CalFormatToImage(triRGBlinCalFormat,Disp.m,Disp.n);

% Converting into XYZ space for Brettel code
load T_xyz1931.mat
T_xyz = SplineCmf(S_xyz1931,T_xyz1931,Disp.wls);
Disp.T_xyz = T_xyz;

% Matrix to convert from rgb to xyz
M_rgb2xyz = Disp.T_xyz*Disp.P_monitor;
M_xyz2rgb = inv(M_rgb2xyz);

triXYZCalFormat = M_rgb2xyz * triRGBlinCalFormat;
triXYZImgFormat = CalFormatToImage(triXYZCalFormat,Disp.m,Disp.n);
[diXYZ] = DichromatSimulateBrettel([], 2, triXYZImgFormat);

diXYZCalFormat = ImageToCalFormat(diXYZ);

diRGBCalFormat = M_xyz2rgb * diXYZCalFormat;
diRGBImgFormat = CalFormatToImage(diRGBCalFormat,Disp.m,Disp.n);
% Quick snipping if only over by a small amount
if max(diRGBImgFormat(:)) < 1.01
    diRGBImgFormat(diRGBImgFormat>1) = .99;
end
diLMSImgFormat = rgbLin2LMSimg(diRGBImgFormat,Disp,1,0);
diLMSCalFormat = ImageToCalFormat(diLMSImgFormat);
% diLMSCalFormat = RGB2LMSimg(diRGBImgFormat,Disp,1,0);

% diLMSCalFormat = simulateDichromatBrettel(lmsImage,renderType,Disp);

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

% CHECK IF MODULATED LMS IS IN GAMUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inGamut = checkGamut(diLMSCalFormat,Disp,bScale);
if inGamut == 0
    error(['tri2dichromatLMSCalFormat: WARNING! rgb values are out of gamut... lmsImage_mod values outside of the range [0 1]']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end