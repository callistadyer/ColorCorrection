function Disp = loadDisplay(img)
% Loads the necessary display parameters.    

if strcmp(img,'gray')
    Disp.m         = 32;
    Disp.n         = 32;
elseif strcmp(img,'ishihara')
    imgSize = 128;
    Disp.m         = imgSize;
    Disp.n         = imgSize;
elseif endsWith(img, '.png', 'IgnoreCase', true) || endsWith(img, '.jpg', 'IgnoreCase', true)
    img_rgb = im2double(imread(img));
    [rows, cols, ~] = size(img_rgb);         
    Disp.m         = cols;
    Disp.n         = rows;
end

% Universal display parameters
wls = (400:10:700)';
d = displayCreate('LCD-Apple');
P_monitor = SplineSrf(displayGet(d, 'wave'), displayGet(d, 'spd'), wls);
load T_cones_ss2;
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,wls);
M_rgb2cones = T_cones*P_monitor;
M_cones2rgb = inv(M_rgb2cones);

Disp.wls       = wls;
Disp.d         = d;
Disp.P_monitor = P_monitor;
Disp.T_cones   = T_cones;

Disp.M_rgb2cones = Disp.T_cones*Disp.P_monitor;
Disp.M_cones2rgb = inv(M_rgb2cones);
Disp.grayRGB     = [0.5 0.5 0.5]';
Disp.grayLMS     = M_rgb2cones*Disp.grayRGB;


end
