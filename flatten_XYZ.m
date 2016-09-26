function [ pointCloudOut ] = flatten_XYZ( XYZ )
%flatten_point_cloud: transforms the input pointcloud by a rotation to a flat
%plane and a translation to the zero plane

%   input: point cloud object
%   output: xyz point cloud flattened according to best fit horizontal
%   plane

% fit planar model to point cloud
[n, ~, ~]           = affine_fit(XYZ);

% create rotation matrix
normalToFlat        = [0 0 1];
rotationAxis        = cross(normalToFlat, n);
rotationAngle       = acos(dot(normalToFlat, n) / ...
                      (norm(normalToFlat) * norm(n))); 
                  
u                   = rotationAxis/norm(rotationAxis);                  
I                   = eye(3);
uCross              = [0, -u(3), u(2); u(3), 0, -u(1),; -u(2), u(1), 0];
uTens               = kron(u,u');    
R                   = cos(rotationAngle)*I + sin(rotationAngle)*uCross ...
                      + (1-cos(rotationAngle))*uTens;
                                  

% rotate the point cloud around the vector defined by the cross product of
% the normal verctor and a vertical vertor by the angle defined by their
% dot product

flatXYZ             = XYZ*R;

% set the data in reference to the average heigth
meanHeight          = nanmean(XYZ(:,3));
zpts                = flatXYZ(:,3)-meanHeight;
xpts                = flatXYZ(:,1)-min(flatXYZ(:,1));
ypts                = flatXYZ(:,2)-min(flatXYZ(:,2));
pointCloudOut       = [xpts,ypts,zpts];

end

