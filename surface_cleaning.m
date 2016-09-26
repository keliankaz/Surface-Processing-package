function [ cleanSurface ] = surface_cleaning(rawSurface, parameterStructure)
%surface_cleaning: converts surface values to nans according to threshold
%parameters
%   Detailed explanation goes here

surface = rawSurface;

if  strcmp(parameterStructure, 'default');
    heightThreshold                 = 8;
    stdThreshold                    = 4;
   
else
    heightThreshold                 = parameterStructure.height;
    stdThreshold                    = parameterStructure.std;
end



%% Height threshold - remove clear outliers is the height field

meanHeight                      = nanmean(nanmean(surface));
heightThreshold                 = meanHeight - ...
                                  heightThreshold*nanstd(nanstd(surface));
tooLow                          = surface < heightThreshold;

outlier                         = tooLow;
surface(outlier)                = nan;

%% model based cleaning
surface                         = fractal_model_outlier(surface,10,stdThreshold);

%% Curvature threshold
surfaceCurvature                = del2(surface);
logSurfaceCurvature             = log(abs(surfaceCurvature));
tooFlat                         = logSurfaceCurvature < -25;
                                  


%% Cleaning
outlier                         = tooFlat;
surface(outlier)                = nan;
cleanSurface                    = surface;
end

