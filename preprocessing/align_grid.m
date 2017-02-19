function newZGrid = align_grid(zGrid, ptSpacing)

rotAngle                = fuckthisshit(zGrid,10000,10);
newZGrid                = rotateZ(zGrid,rotAngle,ptSpacing);
end


function [angle] = fuckthisshit(zGrid, downSampleSize, numScale)

% define downsampling factor in reference to the target downsampled grid
% size.

downSampFactor = ceil(sqrt(numel(zGrid)/downSampleSize));

% downsample 
if downSampFactor >= 1;
    downSampGrid = downsample((downsample(zGrid, downSampFactor))', ...
                          downSampFactor)';
else
    downSampGrid = zGrid;
end

% remove entire rows of nan
[~, M] = size(downSampGrid);
downSampGrid = downSampGrid(sum(isnan(downSampGrid),2)~= M, :);

% define the frequency vector that will be sampled on the rotated grids
[yLength,xLength] = size(downSampGrid);
fVec = logspace(log10(10/(xLength)),log10(1/(2)),numScale);

% unwrap the grid for subsequent rotation

[Xgrid, Ygrid] = meshgrid(1:xLength,1:yLength);

X           = Xgrid(:);
Y           = Ygrid(:);
Z           = downSampGrid(:);
XYZ         = [X,Y,Z];

% remove nan entries in the XYZ points (they prevent regridding)
nanLoc                  = isnan(XYZ(:,3));
XYZ                     = XYZ(~nanLoc,:);

% reggrid at all angles
maxRotation = 180; % degrees
roughness   = zeros(1,maxRotation);

for iAng = 1:maxRotation
    
    rotAngle   = iAng*pi/180;
    rotatedZGrid            = rotateZ(XYZ,rotAngle,1);
    goodRowIndex            = sum(~isnan(rotatedZGrid),2) > 5;
    rotatedZGrid            = rotatedZGrid(goodRowIndex,:); 
    
    % run the periodogram (FFT would probably be better - faster, but this
    % implementation easily takes care of nan entries and avoids alot of
    % hastle.
    
    [PXX, ~] = fault_spectral_density_simple(rotatedZGrid,1,numScale,fVec);
    
    % evaluate the roughness over the entire bandwitdh. Powers are 'log
    % weighted so as not give too large weight to larger scales.
    
    roughness(iAng) = sum(nanmean(log(PXX)));
 
end

% min the amplitude to the peridogram
angle = find(roughness == min(roughness))*pi/180;

end
