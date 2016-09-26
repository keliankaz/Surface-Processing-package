function rotatedGrid = rotateZ(gridOrXYZ,rotAngle,ptSpacing)

% the function assumes that is the column dimension is 3, the data has an
% XYZ form and otherwise it is interpreted as a grid. Rotation follows
% right hand rule convetion around the Z axis.

[yLength,xLength] = size(gridOrXYZ);

if xLength ~=3
    
    [Xgrid, Ygrid] = meshgrid((1:xLength)*ptSpacing,(1:yLength)*ptSpacing);
    
    X           = Xgrid(:);
    Y           = Ygrid(:);
    Z           = gridOrXYZ(:);
    XYZ         = [X,Y,Z];
    
    % remove nan entries in the XYZ points (they prevent regridding)
    nanLoc      = isnan(XYZ(:,3));
    XYZ         = XYZ(~nanLoc,:);
    
else
    XYZ         = gridOrXYZ;
end

% created rotation matrix
Rz                      = [cos(rotAngle), -sin(rotAngle), 0; ...
                           sin(rotAngle),  cos(rotAngle), 0; ...
                           0            ,    0          , 1];
% rotate XYZ points
XYZaligned              = XYZ*Rz;
% redefine x and y grid
newX                    = XYZaligned(:,1);
newY                    = XYZaligned(:,2);
newZ                    = XYZaligned(:,3);

minX                    = min(newX);
minY                    = min(newY);

newX                    = newX-minX;
newY                    = newY-minY;

xVec                    = 0:ptSpacing:max(newX);
yVec                    = 0:ptSpacing:max(newY);
[newXGrid,newYGrid]     = meshgrid(xVec,yVec);

rotatedGrid            = griddata(newX,newY,newZ,newXGrid,newYGrid);

end
