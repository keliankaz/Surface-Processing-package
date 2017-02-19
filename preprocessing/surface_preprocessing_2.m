function [surface, zGrid, ptSpacing ] = surface_preprocessing_2(fileName, unit, instrument)
%surface_preprocessing: takes raw data and converts it into workable
%aligned and gridded output according to desired point (spacing in meters)

%   Detailed explanation goes here

%% load data

disp(fileName)
DATA             = load(fileName);
xyzDATA          = DATA(:,1:3);

% convert into meters

if strcmp(unit, 'mm')
    xyzDATA = xyzDATA/1000;
elseif strcmp(unit, 'cm')
    xyzDATA = xyzDATA/100;
elseif strcmp(unit,'inches')
    xyzDATA = xyzDATA*0.0254;
end

%% remove patches of saturation (for white light) - wont really matter otherwise

if strcmp(instrument,'white light')
% saturated points have a fixed value

saturation      = xyzDATA(:,3) == max(xyzDATA(:,3));
xyzDATA         = xyzDATA(~saturation,:); 

end

 %% flatten data (detrend)

 % remove clear height field oultliers which may latter affect the quality
 % of the fit:
 
 goodData       = xyzDATA(:,3)<(mean(xyzDATA(:,3))+7*std(xyzDATA(:,3)))| ...
                  xyzDATA(:,3)>(mean(xyzDATA(:,3))-7*std(xyzDATA(:,3)));
 xyzDATA        = xyzDATA(goodData,:);
 
 
flatXYZ         = flatten_XYZ(xyzDATA);

%% regrid onto x-y grid

% determine point spacing based on average point density
X               = flatXYZ(:,1);
Y               = flatXYZ(:,2);
Z               = flatXYZ(:,3);

numPts          = length(X);

hullInd         = convhull(X,Y);
hullArea        = polyarea(X(hullInd),Y(hullInd));

ptDensity       = numPts/hullArea;
ptSpacing       = sqrt(1/ptDensity);

% make grid

% create mesh
xAxis           = min(X): ptSpacing : max(X);
yAxis           = min(Y): ptSpacing : max(Y);
[xGrid,yGrid]   = meshgrid(xAxis, yAxis);

surface         = griddata(X,Y,Z,xGrid,yGrid);

%% orient along slickenlines 

% (must come before cleaning because of interpolation steps):

zGrid           = align_grid(surface, ptSpacing);

%% Clean surface
zGrid           = surface_cleaning(zGrid,'default');

end

