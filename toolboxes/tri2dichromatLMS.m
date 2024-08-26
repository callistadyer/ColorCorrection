function [dichromatLMS] = tri2dichromatLMS(lmsImageCalFormat,renderType,cone_mean_orig)

% function takes in trichromat lms values and converts them into dichromat lms values
%
% inputs
% lmsImageCalFormat: lms values of the image for a trichromat
% renderType:        type of dichromacy
%                    Deuteranopia
%                    Protanopia
%                    Tritanopia

l_cone = lmsImageCalFormat(1,:);
m_cone = lmsImageCalFormat(2,:);
s_cone = lmsImageCalFormat(3,:);

% NOTE WELL:
% To get the renderings to look right, we are using the mean value of
% at the missing cone plane to get the scale of our substitution
% right.  Presumably this could be done once by analyzing a full ensemble
% of images, but we are going to worry about that later.

% Cannot rely on 
% deuterMFromLScale = mean(m_cone)/mean(l_cone);
deuterMFromLScale = cone_mean_orig/mean(l_cone);

protoLFromMScale  = mean(l_cone)/mean(m_cone);
tritanSFromMScale = mean(s_cone)/mean(m_cone); 
tritanSFromLScale = mean(s_cone)/mean(l_cone);

% Make dichromat manipulation - missing cone
switch (renderType)
    case 'Deuteranopia' % m cone deficiency
        lmsImageCalFormat(2,:)       =  deuterMFromLScale * l_cone; % replace M cones with L cone PLATE
    case 'Protanopia'   % l cone deficiency
        lmsImageCalFormat(1,:)       = protoLFromMScale*m_cone; % replace L cones with M cone PLATE
    case 'Tritanopia'   % s cone deficiency
        lmsImageCalFormat(3,:)       = (tritanSFromMScale*m_cone + tritanSFromLScale*l_cone)/2; % replace S cones with M cone PLATE
end

dichromatLMS = lmsImageCalFormat;

end