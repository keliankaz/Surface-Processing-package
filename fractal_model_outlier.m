function [cleanSurface] = fractal_model_outlier(surface,numScale, threshold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a model based fractal filter. Functions by rejecting segments with
% unusually large or small RMS values for a given scale 

% Input: 

% surface:      N by M scalar field
% numScale:     Number of scales analysed (scales are logspaced)
% threshold:    Number of standard deviation defining the threshold for
%               defining an outlier.

% Output:

% cleanSurface: The N by M input 'surface' array with outlying segments
%               removed and replaced with nan entires

% The script enssentially loops over every row and collects the RMS,
% location and length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outlierLoc1  = in_one_direction(surface,numScale,threshold);
outlierLoc2  = in_one_direction(surface',numScale,threshold);

surface(outlierLoc1|outlierLoc2') = nan;
cleanSurface = surface;
end

function outlierLoc = in_one_direction(surface,numScale,threshold)

% Defining dimensions
[numy, numx] = size(surface);

% Definition of the scale of observation:
scale0      = 4;                    % smallest scale (cant be less than 3)
maxScale    = numx/2;               % largest scale, numx is num columns
scaleDef    = ceil(logspace(log10(scale0),log10(maxScale),numScale)); 

% initializing locations of oultlying segments
outlierLoc  = zeros(numy,numx);

% working across rows 
for s = 1:numScale % cycle through all the increments of scale
    
    scale       = scaleDef(s);
    
    % initialization to potential max number of entries for segment paramaters:
    maxLength   = numy*(numx-scale+1);
    stdevs      = zeros(1,maxLength);
    stdIndex    = zeros(2,maxLength);
    
    segEntry    = 0;                    % will be useful for storing information in in asym and BFasym
    
% loop over transects parallel to slip
    for ny = 1:numy % number of row
        
        nx = 1;
        
        while nx <= (numx-scale) % nx used to step across the elements in the row segment
                        
            % Defining the segment according to the desired scale
            % Segments are filtered to remove ouliers and stored
            segrange    = (0:scale-1)+nx;           % indexes for the segment
            zsegment    = surface(ny,segrange);     % the topography transect segment
            
            % checking for entries that are 'nan'
            znancheck   = isnan(zsegment);
            
            if sum(znancheck)==0% making sure that there are no nan in the segment
                 
                 % advacing counter
                 segEntry   = segEntry + 1; % number of looped segments
                 
                 % Saving Standard deviations on segment:
                 stdevs(1,segEntry) = std(detrend(zsegment));
                 stdIndex(:,segEntry) = [nx;ny];
                 
                 %advancing counter
                 nx = nx+1;
            else
                 % Skip all the next segments that will include a nan:
                 ind        = find(znancheck,1,'last');
                 nx         = nx + ind;
            end        
        end
    end
    
    % This conditional statement stops the code when there are no longer
    % segments that are long enough to study
    if segEntry == 0
        break
    else 
        % Extracting only non-sparse elements       
        stdevs               = stdevs(1,1:segEntry);
        
        % Identify Location of unwanted Segments
        outlierSegmentLoc    = stdIndex(:, ...
                              (stdevs < mean(stdevs) - threshold*std(stdevs) | ...
                               stdevs > mean(stdevs) + threshold*std(stdevs)));
                                    
        % identify entries that will be removed by turning bad locations
        % from 1 to 0 by enry-wise multiplication. 
        
        for iOut = 1:length(outlierSegmentLoc)
            xLoc                = outlierSegmentLoc(1,iOut);
            yLoc                = outlierSegmentLoc(2,iOut);
            segrange            = (0:scale-1) + xLoc;
            outlierLoc(yLoc,segrange) = 1;
        end
    end
end
end