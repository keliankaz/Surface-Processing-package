function [] = unpack_parameters(desiredPlot, varargin)
% unpacks the output from the surface processing code package.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT:

% desiredPlot:
%   'FFT'           : power spectrum as determined by the fast fourier analysis
%   'PLOMB'         : periodogram as determined by the plomb analysis
%   'avgAsyms'      : scale dependent average asymetry
%   'topostd'       : RMS plot
%   'topoSkew'      : scale dependent skewness of height fields
%   'topoKurt'      : scale dependent kurtoisis of height fields
%   'PowerVsDisp'   : power at a given scale as a function of displacement
%   'RMSVsDisp'     : model RMS at a given scale as a function of
%                     displacement
%   'Grids'         : shows both the original and pre-processed grid for
%                     the specified file 'fileName'
%   'Best Fits'     : best logarythmic fits to power spectra obtained from
%                     FFT

% varargin:
%   orientation:
%       'parallel'      : slip parallel analysis
%       'perpendicular  : slip perpendicular analysis
%   displacement    :(optional) array with desplacements in the corresponding
%                     order to the input files.
%   scale           : required for 'powerVsDisp'. Also requires displcament
%                     array to be included. Specifies the scales at which
%                     scale the 'powerVsDisp is going to be plotted.
%                     Interpolation is done accoding to the linear
%                     regression model for the entire PLOMB spectrum
%   Constraint      : Cell array specifying constraint on displacement 
%                     ('Upper Bound' or 'Direct')
%   parameter (with 'PowerVsDisp',orientation, displacement, scale, constraint, ...):
%       'Hurst'         : see the evolution of the Husrt exponent with
%                         displacement
%       'prefactor'     : evolution of the prefactor with displacement
%   fileName        : string with the desired file name for the 'Grids' plot



% User will be prompted to navigate to the directory in which the surface
% processing output files are stored. The folder must not have anyother
% files inside it.

% OUTPUT:
%   plot of the specified parameters as a function of scale (or displacement)
%   for all input files in directory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% roughness at a given scale as a function of displacement

inputs = varargin;

if  strcmp(desiredPlot,'PowerVsDisp') || strcmp(desiredPlot,'RMSVsDisp')
    roughnessVsDisp(desiredPlot, inputs) 
else 
    if strcmp(desiredPlot, 'Grids')
        plotgrid(inputs)
    else
        plotspectra(desiredPlot, inputs)
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% save plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% location = 'C:\Users\kelian_dc\Documents\school\Masters thesis\Thesis\Data_processing\Figures\New PLots'; % desired save location
% saveFileName = [location,'\','LiDar',' - ',desiredPlot,'-',inputs{1},'-',date];
% savefig(saveFileName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


%% plot a metric of roughness vs displacement
function [] = roughnessVsDisp(desiredPlot, inputs)

[files, numFiles,fileIndex] = getfiles();

% identify data
orientation  = inputs{1};
displacement = inputs{2};
scale        = inputs{3};

numIn = length(inputs);

if numIn >= 4
    constraint = inputs{4};
end

if numIn >= 5
    wellProcessed = inputs{5};
end

if numIn >= 6
    parameter = inputs{6};
end

% create the displacement "domain" - for line plotting

minD        = min(displacement(displacement~=0));
maxD        = max(displacement);
minmaxD     = [minD,maxD];

Power       = zeros(1,numFiles);
colorSpec   = zeros(numFiles,3);
fileNameArray = cell(1,numFiles);

for iFile = 1:numFiles
    
    fileName            = files(fileIndex(iFile)).name;
    load(fileName,'parameters')
    fileNameArray{1,iFile} = getfield(parameters,'fileName');
    spectrumType        = 'FFT';
    desiredData         = getfield(getfield(parameters, orientation),spectrumType);
    fx                  = desiredData{1,1};
    Px                  = desiredData{1,2};
    Px                  = Px(1:length(fx));
    PxNanInd            = isnan(Px);
    fx                  = fx(~PxNanInd);
    Px                  = Px(~PxNanInd);
    
    d                   = makebestfit(fx,Px);
%     d                   = polyfit(log10(1./fx),log10(Px'),1);
    
    % this is a bit awkward but fuckit im tired (if nargin is 7 we
    % basically continue with the same code but instead of using the
    % power we just use the Hurst exponent (the slope of the best fit
    % line through the data).
    
    if numIn == 6
        if strcmp(parameter, 'Hurst')
            Power(iFile)        = (d(1)-1)/2;
        else strcmp(parameter, 'prefactor')
            Power(iFile)        = d(2);
        end
    elseif strcmp(desiredPlot,'PowerVsDisp')
        Power(iFile)        = d(1)*log10(scale) + d(2);
        
    % Handle the conversion to RMS as done in Brodsky et al., 2011
    % P(lambda) = C*lambda^BETA
     
    elseif strcmp(desiredPlot,'RMSVsDisp')
        C = 10^d(2);
        BETA = d(1);
        RMS = (C/(BETA - 1))^0.5*scale*(BETA-1)/2;
        Power(iFile) = log10(RMS); %... for the purpose of plotting...
    end
    
    % identify constraints on data points
    if strcmp(constraint{iFile},'Upper Bound')
        colorSpec(iFile,:) = [0.5 0.5 0.5];
    elseif strcmp(constraint{iFile}, 'Direct')
        colorSpec(iFile,:) = [0 0 0];
    end
end

zeroInd             = displacement == 0;
zeroDispPower       = Power(zeroInd);

figure
scatter(log10(displacement),Power,40,colorSpec, 'filled')
hold on
xlabel('log(displacement)')
if numIn == 6
    ylabel(parameter)
    title([parameter, ' as a function of displacement at ',num2str(scale),' using ', spectrumType,' - ' date])
else
    ylabel('log(Power)')
    title(['Roughness as a function of displacement at ',num2str(scale),' using ', spectrumType,' - ' date])
end

% tags to the points for analysis

if numIn == 5
    noLoc = strcmp(wellProcessed, 'no');
    colorSpec = zeros(numFiles,3);
    for iFile = 1:numFiles
        if noLoc(iFile) == 1
            colorSpec(iFile,:) = [1 0 0];
        end
    end
else
    colorSpec = zeros(numFiles,3);
end

% offset = 0.01;
% % t = text(log10(displacement)+offset,Power+offset,fileNameArray,'interpreter', 'none');
% 
% for iFile = 1:numFiles
%     t(iFile).Color = colorSpec(iFile,:);
% end



% Zero displacement Data ploted as lines
for iLine = 1:length(zeroDispPower)
    plot(log10(minmaxD), zeroDispPower(iLine)*[1,1], 'r')
end

nanInd    = isnan(displacement);
goodInd   = ~nanInd & ~zeroInd;
goodData  = displacement(goodInd)';
goodPower = Power(goodInd);

% pass a best fit line through the entire dataset
d = polyfit(log10(goodData),goodPower,1);
plot(log10(minmaxD),(d(1)*log10(minmaxD) + d(2)));

% pass a best fit through well constrained data
constrainedInd = strcmp(constraint,'Direct');
goodInd   = constrainedInd & ~zeroInd;
goodData  = displacement(goodInd)';
goodPower = Power(goodInd);

% pass a best fit line through the entire dataset
d = polyfit(log10(goodData),goodPower,1);
plot(log10(minmaxD),(d(1)*log10(minmaxD) + d(2)), 'Linewidth',2);

hold off
end

%% plot a specific, or many, grids (surfaces) 
function [] = plotgrid(inputs)

desiredFileNameArray = inputs{1};
numGrid  = length(desiredFileNameArray);
subplotCount = 0;

for iGrid = 1:numGrid;
    
    load(desiredFileNameArray{iGrid})
    originalGrid = getfield(parameters,'zGrid');
    cleanGrid    = getfield(parameters,'newZGrid');
    
    subplotCount = subplotCount + 1;
    
    subplot(iGrid,2,subplotCount)
    imagesc(originalGrid)
    axis equal
    xlabel('x')
    ylabel('y')
    title('Original Grid')
    
    subplotCount = subplotCount + 1;
    
    subplot(numGrid,2,subplotCount)
    imagesc(cleanGrid)
    axis equal
    xlabel('Slip Perpendicular')
    ylabel('Slip Parallel')
    titel('Pre-processed Grid')
end
end

%% plot all the frequency spectra in a given file
function [] = plotspectra(desiredPlot,inputs)

[files, numFiles,fileIndex] = getfiles();

legendArray = cell(1,numFiles);
orientation = inputs{1};
hold on

if length(inputs) >= 2
    displacement = inputs{2};
    D = displacement;
    minD = min(D(D~=0));
    maxD = max(D);
    logMinD = log10(minD);
    logMaxD = log10(maxD);
    logDelD = logMaxD-logMinD;
    
end

if length(inputs) >=3
    constraint = inputs{3};    
end

% loop over files
for iFile = 1:numFiles
    fileName            = files(fileIndex(iFile)).name;
    load(fileName)
    
    if strcmp(desiredPlot,'best fits')
        desiredStruct   = getfield(getfield(parameters, orientation),'FFT');
        desiredCell     = desiredStruct;
        
        wavelength      = 1./desiredCell{1};
        power           = desiredCell{2};
        
        power           = power(1:length(wavelength));
        nanInd          = isnan(power);
        power           = power(~nanInd);
        wavelength      = wavelength(~nanInd);
        
        d               = polyfit(log10(wavelength'),log10(power),1);
        x               = [min(wavelength),max(wavelength)];
        y               = 10.^(d(1)*log10(x)+d(2));
        
    else
        desiredStruct   = getfield(getfield(parameters, orientation),desiredPlot);
        
        % decide what plot to do (enable plot-specific attributes)
        if strcmp(desiredPlot,'FFT') || strcmp(desiredPlot,'PLOMB')
            desiredCell     = desiredStruct;
            x               = desiredCell{1};
            x               = 1./x;
            y               = desiredCell{2};
            y               = y(1:length(x));
            
        else
            x               = getfield(getfield(parameters, orientation),'Scales');
            if strcmp(desiredPlot, 'avgAsyms')
                badInd          = isnan(desiredStruct)|desiredStruct == 0;
                y               = abs(desiredStruct(~badInd));
                x               = x(~badInd);
            else
                y               = desiredStruct;
                y               = y(y~=0);
                x               = x(y~=0);
            end
            
        end
    end
    
    scanName = getfield(parameters,'fileName');
   
    % create the plots
    if length(inputs) == 1
        plot(x,y,'.-')
        legendArray{1,iFile} = fileName;
        
    elseif length(inputs) >= 2
        % set color as a function of displacement
        if D(iFile) == 0
            colorForDisp = [0 1 0];
        else
            colorInd = (log10(D(iFile))-logMinD)/(2*(logDelD));
            colorForDisp = [0.5-colorInd, 0.5-colorInd, 0.5-colorInd];
        end
        
        % set unknow displacements to grey
        if sum(isnan(colorForDisp)) ~= 0
            colorForDisp = [0.5, 0.5, 0.5];
        end
        
        p = plot(x,y,'-');
        p.Color = colorForDisp;
        formatSpec = 'D = %d - file: %s';
        legendArray{1,iFile} = sprintf(formatSpec, D(iFile), ...
            scanName);
        xlabel('Scale (m)')
    end
    
     % to help distinguish points whith various constraints on displacement
    if length(inputs) >=3
        if strcmp(constraint{iFile}, 'Direct')
            p.LineWidth = 3;
        elseif strcmp(constraint{iFile},'Upper Bound')
            p.LineWidth = 0.2;
        end
    end
   
    
    ylabel(desiredPlot)
    title([desiredPlot,' as a function of scale - ', orientation, date])
    
    if strcmp(desiredPlot,'topostd') || strcmp(desiredPlot,'FFT')|| ... 
            strcmp(desiredPlot,'PLOMB')
        set(gca,'XScale','log','YScale','log')
        xlabel('Scale (log(m))')
    end  
end

legend(legendArray,'interpreter', 'none')

end

%% small function to query the files in the desired directory
function [files, numFiles, fileIndex] = getfiles()
directory_name = uigetdir;
files = dir(directory_name);
addpath(directory_name)

fileIndex = find(~[files.isdir]& ~ strcmp({files.name},'.DS_Store'));
numFiles  = length(fileIndex);
end

%% best fit the spectra while not letting instrument artefact affect results
function [d] = makebestfit(F,P,varargin);
    % small function to make the best fit of the power spectrum
    
    % ensure the vectors are of the same size/dimension
    if length(F) == length(P)
        if size(F) ~= size(P)
            P = P';
        end
        
        % make default setting
        defaultStruct.FitMethod     = 'default';
        defaultStruct.SectionVal    = 1/2;
        
        % adapt defualt structure
        setVal('FitMethod')
        setVal('SectionVal')
        
        S = defaultStruct;

        if strcmp(S.FitMethod, 'default')
            % Simplest form: take the fit over the entire calculated spectrum:
            d           = polyfit(log10(1./F),log10(P),1);
        
        elseif strcmp(S.FitMethod, 'section')
                    % Uninspired, but better form: take the best fit only on a
        % specifict section of the fit.
        
        % In the case of the scan data, this
        % only takes into account the larger wave lengths where the data
        % should not be affected by instrumental artefacts
            
            numIn       = length(F);
            numEntries  = ceil(numin*S.SectionVal);
            start       = numIn-numEntries;
            newF        = F(start:end);
            newP        = P(start:end);
            
            d           = polyfit(log10(1./newF,newP));
        end
        
    
    else
        disp('error in the code making the best fit, input vectors must be the same length')
    end
    
    
    function setVal(name)
        ind = strcmp(name, varargin);
        if any(ind)
            setfield(defaultStruct,name,varargin{1,ind+1});
        end
    end
        
end

    
    
