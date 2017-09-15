function [FITCOEFF,FITERROR] = runmontecalrofit(N,DATA, ERROR, varargin)
% run motecarlos simulation to obtain best fit parameters based on all data
% points, errors both in x and y 

% input:
% N:        number of iterations
% DATA:     data going into the fit (e.g. [X,Y]
% ERROR:    corresponding 1 sigma error (according to error model)

% varargin: (input pairs)
% ...,xErrorModel,ERRORMODEL,... : type of error model 'gaussian', 'lognormal'
%                                  or 'equal'. Cell array with dimension
%                                  equivalent to length(DATA)
% ...,fitModel,FITMODEL,...      : out of use do not use
% ...,histogram, 'on',...        : plot histogram of results 'on' or 'off'
%                                 (default)

% not ideal here:  remove rows that contain nan
ind = logical(sum(DATA~=DATA,2));
DATA(ind,:)     =[];
ERROR(ind,:)    =[];

dim             = size(DATA);
numPoints       = max(dim); % number of points in the data set
numDim          = min(dim); % number of dimensions in the dataset

userInput       = {'errorModel'     , ...
                   'fitModel'       , ...
                   'histogram'      };
                     
defaultStruct               = [];
defaultStruct.errorModel    = repmat({'gaussian'},numDim,numPoints);
defaultStruct.histogram     = 'off';

defaultStruct.fitModel      = 'power1';

S = setVal(defaultStruct,userInput, varargin);

% initialize values:
errorModelTypes = unique(S.errorModel(:));

fitType     = fittype(S.fitModel);
numCoeffs     = numcoeffs(fitType);

fitCoeffs   = zeros(N,numCoeffs);


% consider parallel processing
for n = 1:N
    
    newDATA     = zeros(numPoints,numDim);
    
    % resample points from their respective distributions
    if      any(strcmp(errorModelTypes,'gaussian')) 
        negVal = 1;
        while any(negVal)
            loc = find(strcmp(S.errorModel,'gaussian'));
            newDATA(loc)    = normrnd(DATA(loc),ERROR(loc));
            negVal = newDATA(loc)<0;
        end
    end 
        
    if  any(strcmp(errorModelTypes,'lognormal'))
        
        loc = find(strcmp(S.errorModel,'lognormal'));
        newDATA(loc)    = lognrnd(log(DATA(loc)),ERROR(loc));
        
    end
    
    if  any(strcmp(errorModelTypes,'equal'))
        
        loc             = find(strcmp(S.errorModel,'equal'));
        numInd          = length(loc);
        newDATA(loc)    = DATA(loc) + (2*rand(numInd,1)-1).*ERROR(loc);
    end
    
    d               = polyfit(log10(newDATA(:,1)),log10(newDATA(:,2)),1);
    fitCoeffs(n,:)  = d;
    
end

FITCOEFF = nanmean(fitCoeffs); % the nan here is not ideal
FITERROR = nanstd(fitCoeffs);  % because it effectively removes the pts with H<1

if strcmp(S.histogram,'on')
    figure
    for iCoeff = 1:numCoeffs
        subplot(1,numCoeffs,iCoeff)
        histogram(fitCoeffs(:,iCoeff))
        xlabel('value')
        ylabel('count')
    end
    disp(fitType)
    disp(FITCOEFF)
    disp(FITERROR)
end
end