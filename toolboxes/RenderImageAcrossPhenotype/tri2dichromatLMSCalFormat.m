function [dichromatLMSCalFormat] = tri2dichromatLMSCalFormat(lmsImageCalFormat,renderType,cone_mean_orig,Disp,bScale)

% function takes in trichromat lms values and converts them into dichromat lms values
%
% Syntax:
%   [dichromatLMS] = tri2dichromatLMSCalFormat(lmsImageCalFormat,renderType,cone_mean_orig,Disp,bScale)
%
% Description:
%
% Inputs:
%   lmsImageCalFormat  LMS values of the image for a trichromat
%   renderType         Type of dichromacy
%                        Deuteranopia
%                        Protanopia
%                        Tritanopia
%   cone_mean_orig     Mean of missing cone values for trichromat (original image)
%                      This is used to calculate missing cone value for dichromat 
%
% Outputs:
%   dichromatLMS       LMS values for dichromat type specified by renderType
%
% Optional key/value pairs:
%   None

l_cone = lmsImageCalFormat(1,:);
m_cone = lmsImageCalFormat(2,:);
s_cone = lmsImageCalFormat(3,:);

% NOTE WELL:
% To get the renderings to look right, we are using the mean value of
% at the missing cone plane to get the scale of our substitution
% right.  Presumably this could be done once by analyzing a full ensemble
% of images, but we are going to worry about that later.

% NOTICE!! THINK THIS IS MESSING WITH THE GAMUT
deuterMFromLScale = cone_mean_orig/mean(l_cone);
protoLFromMScale  = cone_mean_orig/mean(m_cone);
tritanSFromMScale = cone_mean_orig/mean(m_cone); 
tritanSFromLScale = cone_mean_orig/mean(l_cone);

% Make dichromat manipulation - missing cone
switch (renderType)
    case 'Deuteranopia' % m cone deficiency
        lmsImageCalFormat(2,:)       = l_cone; 
        [scaleFactor,dichromatLMSCalFormat,inGamut] = findDichromMappingScaleFactor(2,lmsImageCalFormat,Disp,bScale);

    case 'Protanopia'   % l cone deficiency
        lmsImageCalFormat(2,:)       = m_cone; 
        [scaleFactor,dichromatLMSCalFormat,inGamut] = findDichromMappingScaleFactor(1,lmsImageCalFormat,Disp,bScale); 
        
    case 'Tritanopia'   % s cone deficiency
        lmsImageCalFormat(3,:)       = (m_cone + l_cone)/2;
        [scaleFactor,dichromatLMSCalFormat,inGamut] = findDichromMappingScaleFactor(3,lmsImageCalFormat,Disp,bScale);

end

%dichromatLMSCalFormat = lmsImageCalFormat;

% CHECK IF MODULATED LMS IS IN GAMUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inGamut = checkGamut(dichromatLMSCalFormat,Disp,bScale);
if inGamut == 0
    error(['tri2dichromatLMSCalFormat: WARNING! rgb values are out of gamut... lmsImage_mod values outside of the range [0 1]']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end