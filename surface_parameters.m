function [parameterStruct] = surface_parameters(surface ,pointSpacing, numberOfScales)
% This function analysises surface roughness metrics: standard deviation,
% skewness, and kurtosis. Specifically it also measures directional
% asymmetries in the surfacea. This is a side project of Kelian
% Dascher-Cousineau's master's thesis. Note that this code is is expensive
% (sensitive to the number of scales specified.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT: 

% pointSpacing is the gridspacing in meters or the surface scan

% numberOfScales is the desired number of sampled scales at which the
% measurements for asymmetry are cummulated (note that scales are then
% logspaced)

% OUTPUT:

% parameterStruct:

%   parameterStruct.Scales      : Sampled scales for parameterization
%   parameterStruct.avgAsyms    : Estimation of scale dependent asymetry
%   parameterStruct.topostd     : RMS as a function of scale
%   parameterDtruct.topoSkew    : skewness in the height field with scale
%   parameterStruct.topoKurt    : kurtosis in the height field with scale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the derivative of the surface

[zp,~]          = gradient(surface,pointSpacing); 
h               = pointSpacing; % point spacing in m

% Defining dimensions
[numy, numx] = size(surface);

% Definition of the scale of observation:

scale0      = 4;                    % smallest scale (cant be less than 3)
maxScale    = numx;                 % largest scale, numx is num columns
numScale    = numberOfScales;       % number of scale increments to include

% defining the vector with all the scales of interest
scaleDef    = ceil(logspace(log10(scale0),log10(maxScale),numScale)); 

% initialising scale dependent average asymetry array:
avgAsyms    = zeros(1,numScale);    
topostd     = zeros(1,numScale);    
topoSkew    = zeros(1,numScale);
topoKurt    = zeros(1,numScale);
        

% initializing cell arrays saving info about each segment for given scales

asym        = cell(numScale,1);
stdevs      = cell(numScale,1);
skewArray   = cell(numScale,1);
kurtArray   = cell(numScale,1);

% working across rows

for s = 1:numScale % cycle through all the increments of scale
    
    scale           = scaleDef(s);
    
    % initialization to potential max number of entries for segment paramaters:
    maxlength       = numy*(numx-scale+1);
    
    asym{s,1}       = zeros(1,maxlength);   % asymetry on individual segments
    stdevs{s,1}     = zeros(1,maxlength);   % stadard deviations on segments
    skewArray{s,1}  = zeros(1,maxlength); 
    kurtArray{s,1}  = zeros(1,maxlength); 
    
    segEntry        = 0;                    % will be useful for storing information in in asym and BFasym
    
% loop over transects parallel to slip
    for ny = 1:numy % number of rows
        
        nx = 1;
        
        while nx <= (numx-scale) % nx used to step across the elements in the row segment
                        
            % Defining the segment according to the desired scale
            % Segments are filtered to remove ouliers and stored
            segrange    = (0:scale-1)+nx;                                               % indexes for the segment
            
            zsegment    = surface(ny,segrange);                                         % the topography transect segment
            zsegment((zsegment > ((nanmean(zsegment))+4*nanstd(zsegment)))) = nan;      % filter noisy parts
            
            zpsegment   = zp(ny,segrange);                                              % the derivative slope transect
            
            % checking for entries that are 'nan'
            znancheck   = isnan(zsegment);        % ... in the topography data
            zpnancheck  = isnan(zpsegment);      % ... in the slope data
            
            if sum(znancheck)==0 && sum(zpnancheck)==0 % making sure that there are no nan in the segment
                 
                 % advacing counter
                 segEntry   = segEntry + 1; % number of looped segments
                
                 % Detrending the segment with respect to least squares
                 % linear slope on the segment:
                 
                 % slope of the segment:
                 xs         = h*segrange;
                 xbar       = mean(xs);
                 zbar       = mean(zsegment);
                 slope      = (sum(xs.*zsegment)-scale*xbar*zbar)/(sum(xs.^2)-scale*(xbar^2));
                 
                 % topographic slope segment is removed from derivative
                 fzp        = zpsegment-slope;
                 detrendedZSeg = detrend(zsegment);
                 
                 % Positive and negative slopes are identified and
                 % averaged:
                 zpPos      = fzp(fzp>0);  
                 zpNeg      = fzp(fzp<0);                   
                 avgzpPos   = mean(zpPos);                
                 avgzpNeg   = mean(zpNeg);                  
                 
                 % assymetry measured as the sum of the average positive
                 % and average negative derivative
                 asym{s,1}(1,segEntry) = avgzpPos + avgzpNeg;
                 
                 % Saving Standard deviations on segment:
                 stdevs{s,1}(1,segEntry) = std(detrendedZSeg);
                 
                 % Saving skewness of height distribution on segment:
                 skewArray{s,1}(1,segEntry) = skewness(detrendedZSeg);
                 
                 % Saving kurtosis of height distribution on segment:
                 kurtArray{s,1}(1,segEntry) = kurtosis(detrendedZSeg);
                 
                 %advancing counter
                 nx = nx+1;
            else
                 % Skip all the next segments that will include a nan:
                 nanind     = find(znancheck==1);
                 nanindp    = find(zpnancheck==1);
                 ind        = max([nanind,nanindp]);
                 nx         = nx + ind;
                 
            end        
        end
    end
    
    % terminate script if there are no continuous segments for the scales
    % lager than the active one
    if segEntry == 0
        break
    else
        % Extracting only non-sparse elements
        asym{s,1}           = asym{s,1}(1,1:segEntry);
        stdevs{s,1}         = stdevs{s,1}(1,1:segEntry);
        skewArray{s,1}      = skewArray{s,1}(1,1:segEntry);
        kurtArray{s,1}      = kurtArray{s,1}(1,1:segEntry);
        
        % Asymmetry recorded for this given scale:
        avgAsyms(s)         = mean(asym{s,1});
        topostd(s)          = mean(stdevs{s,1});
        topoSkew(s)         = mean(skewArray{s,1});
        topoKurt(s)         = mean(kurtArray{s,1});
        
    end
end

% Storing surface parameters into a structure
parameterStruct.Scales      = scaleDef*h;
parameterStruct.avgAsyms    = avgAsyms;
parameterStruct.topostd     = topostd;
parameterStruct.topoSkew    = topoSkew;
parameterStruct.topoKurt    = topoKurt;

end

