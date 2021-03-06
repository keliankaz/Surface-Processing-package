function [grid, pointSpacing] = parse_zygo_format(varargin)
% pasrse through the zygo export format and create a 2D grid of the
% topographic data and also outputs the point spacing

% input must be in pair: 'fileName' followed by the file name (string) -
% required input field 
% 'detrend' followed by 'yes' or 'no'(default is 'no') - optional field
% (NOTE this is detrending NOT rotating - DO NOT USE FOR VERY INCLINED
% SURFACES!!!!)

% possible vargargin
userInput   = {'fileName', ...
               'detrend',...
               'infill'};
               
defaultStruct           = [];
defaultStruct.detrend   = 'no';
defaultStruct.infill    = 'yes';

S = setVal(defaultStruct, userInput, varargin);

delimeter = ' ';
S.fileName
fileName = S.fileName;
FID       = fopen(fileName);

numY = dlmread(fileName,delimeter,[3 2 3 2]);
numX = dlmread(fileName,delimeter,[3 3 3 3]);
pointSpacing = dlmread(fileName,delimeter,[7 6 7 6]);

Zdata = textscan(FID,'%*f %*f %f','Delimiter',' ','HeaderLines', ...
                 14,'TreatAsEmpty','No Data');
z     = Zdata{1,1};
z     = z*10^-6; % convert to meters

grid = reshape(z,numY,numX);

if strcmp(S.detrend,'yes')
    
    % creat x-y coordinates
    xVec = 1:numX;
    yVec = 1:numY;
    [X,Y] = meshgrid(xVec,yVec);
    
    x = X(:);
    y = Y(:);
    
    % remove nan values
    goodInd = ~isnan(z);
    x       = x(goodInd);
    y       = y(goodInd);
    z       = z(goodInd); 
        
    meanZ   = mean(z);
    z       = z-meanZ;
  
    [n,~,p] = affine_fit([x,y,z]); % compute best fitting plane equation
    
    % remove the mean from the grid
    grid    = grid - meanZ; 
    
    % remove planar trend (NOTE this is detrending NOT rotating - DO NOT
    % USE FOR VERY INCLINED SURFACES!!!!)
    
    planeTrend = -(n(1)*(X-p(1)) + n(2)*(Y-p(2)))/n(3)+p(3); 
    grid = grid - planeTrend;
end

if strcmp(S.infill,'yes')
   %  use a the inpaint nan inteprolation scheme with the spring metaphor
   %  (method 4). Function from the file exchange:
   %  https://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint-nans
   
   grid     = inpaint_nans(grid);
   
end
    

end