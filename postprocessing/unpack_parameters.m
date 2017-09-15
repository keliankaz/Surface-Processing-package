function [] = unpack_parameters(desiredData, varargin)
% unpacks the output from the surface processing code package.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT:

% desiredData:
%   'FFT'           : power spectrum as determined by the fast fourier analysis
%   'PLOMB'         : periodogram as determined by the plomb analysis
%   'avgAsyms'      : scale dependent average asymetry
%   'topostd'       : RMS plot
%   'topoSkew'      : scale dependent skewness of height fields
%   'topoKurt'      : scale dependent kurtoisis of height fields
%   'PowerVsDisp'   : power at a given scale as a function of displacement
%   'RMSVsDisp'     : model RMS at a given scale as a function of
%                     displacement
%   'DirectRMSVsDisp: direct measurements of RMS as a function of
%                     displacement.
%   'Grids'         : shows both the original and pre-processed grid for
%                     the specified file 'fileName'
%   'Best Fits'     : best logarythmic fits to power spectra obtained from
%                     FFT
%   'hurstVsPrefactor': correlation between best fit hurst exponent and
%                     prefactor
%   

% varargin:
%   orientation:
%       'parallel'      : slip parallel analysis
%       'perpendicular  : slip perpendicular analysis
%   scale           : required for 'powerVsDisp'. Also requires displcament
%                     array to be included. Specifies the scales at which
%                     scale the 'powerVsDisp is going to be plotted.
%                     Interpolation is done accoding to the linear
%                     regression model for the entire PLOMB spectrum
%   displacement    :(optional) array with desplacements in the corresponding
%                     order to the input files.
%   Constraint      : Cell array specifying constraint on displacement 
%                     ('Upper Bound' or 'Direct')
%   parameter (with 'PowerVsDisp',orientation, displacement, scale, constraint, ...):
%       'Hurst'         : see the evolution of the Husrt exponent with
%                         displacement
%       'prefactor'     : evolution of the prefactor with displacement
%   fileName        : string with the desired file name for the 'Grids' plot

%   pair-wise inputs:
%       'magnification',mag: specifier magnification followed by the
%                            desired magnification
%       'occular', occular:  "                                                     
%                                                 "
%       'subset',{specifier,subsets}: choose only a subset of all scans
%                                     based on a field characteristic or
%                                     value (e.g. ...,'subset',{magnification,
%                                     [10,20]},...)


% User will be prompted to navigate to the directory in which the surface
% processing output files are stored. The folder must not have anyother
% files inside it.

% OUTPUT:
%   plot of the specified parameters as a function of scale (or displacement)
%   for all input files in directory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% examples:

%% 1.
% To plot the Power exponent as a function of displacement ('PowerVsDisp'),
% at a scale of 0.0001 ('scale' ,0.0001), only considering scans taken with
% the 2.75 occular ('scale' ,0.0001), with sample names displayed next to
% points for reference ('text', 'on'):

% unpack_parameters('roughnessVsDisp','parameter','Power','scale' ,0.0001,'subset',{'occular',2.75},'text', 'on')

%% 2. 
% To plot all the fast fourier transform ('FFT') spectra in a directory
% color coded according to the displacement values specified in their
% respective -parameter- structure, in the parallel direction
% ('orientation','parallel'):

% unpack_parameters('FFT','orientation','parallel')

%% 3.
% To plot all RMS(L) curves in a directory color coded according to the
% displacement values specified in their respective -parameter- structure,
% in the parallel direction ('orientation','parallel'):

% unpack_parameters('topostd','orientation','parallel')

%% 4.
% Print a list of a given parameter (Husrt exponent, prefactor, power or
% RMS at a given scale simply:

% unpack_parameters('Hurst')
% unpack_parameters('Power','scale', 0.25)
% unpack_paramaters('RMS_internet_stylez_with_a_Z', 'scale', 0.0002) % compute RMS by integrating under the FFT sprectrum

%% feed input into correct workflow (functions)

if ~strcmp(desiredData,'Grids')
    [files, numFiles,fileIndex] = getfiles();
    [S,subsetLoc,numFiles]= parseInput(varargin,files,numFiles,fileIndex); 
end

if          strcmp(desiredData, 'roughnessVsDisp'    ); roughnessVsDisp         (S,files, fileIndex, subsetLoc,numFiles); 
elseif      strcmp(desiredData, 'Grids'              ); plotgrid(varargin);
elseif      strcmp(desiredData, 'hurstVsPrefactor'   ); hurstVsPrefactor        (S,files, fileIndex, subsetLoc,numFiles);
elseif  any(strcmp(desiredData,{'Hurst'     , ...
                                'prefactor' , ...
                                'Power'     , ...
                                'RMS_emily' , ...
                                'RMS_internet_stylez_with_a_Z' ...
                                                   })); get_suface_parameter    (desiredData, S, files, fileIndex ,subsetLoc,numFiles,'yes');     
else                                                  ; plotspectra(desiredData, S, files, fileIndex ,subsetLoc,numFiles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% save plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% location = 'C:\Users\kelian_dc\Documents\school\Masters thesis\Thesis\Data_processing\Figures\New PLots'; % desired save location
% saveFileName = [location,'\','LiDar',' - ',desiredPlot,'-',inputs{1},'-',date];
% savefig(saveFileName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%% Main function
%% plot a metric of roughness vs displacement
    function [] = roughnessVsDisp(S, files, fileIndex, subsetLoc,numFiles)
    
% extract all the info from the strucutre now that it is set 
parameter           = S.parameter;     
scale               = S.scale;   
displacement        = S.displacement;
constraint          = S.constraint;
wellProcessed       = S.wellProcessed;
displacementError   = S.displacementError;

% chack if displacement is defined for any of the scans, in none are
% defined, this function makes no sense

if isempty(displacement)
    error('Error using roughnessVsDis, displacement values must be specified on at least one of the scans');
end

% create plots:

% create the displacement "domain" - for line plotting

minD        = min(displacement(displacement~=0));
maxD        = max(displacement);
minmaxD     = [minD,maxD];

[fileNameArray, Power, powerError] = get_suface_parameter(parameter, S, files, fileIndex, subsetLoc,numFiles,'no');

% clean out data
nanInd    = isnan(displacement);
zeroInd   = displacement == 0;

% make data same dimension
if size(displacement) ~= size(Power); Power = Power'; end

% classify data

allConstraint   = string(constraint        (~zeroInd & ~nanInd)');
allDisp         = displacement      (~zeroInd & ~nanInd);
allDispError    = displacementError (~zeroInd & ~nanInd);
allPower        = Power             (~zeroInd & ~nanInd);
allPowerError   = powerError        (~zeroInd & ~nanInd, :);

zeroDispPower       = Power(zeroInd);

upperBoundInd       = strcmp(allConstraint,'Upper Bound');
upperBoundPower     = allPower         (upperBoundInd);
upperBoundPowerErr  = allPowerError    (upperBoundInd, : );
upperBoundDisp      = allDisp          (upperBoundInd);
upperBoundDispError = allDispError     (upperBoundInd);

directInd           = strcmp(allConstraint,'Direct');
directPower         = allPower         (directInd);
directPowerErr      = allPowerError    (directInd, :);
directDisp          = allDisp          (directInd); 
directDispError     = allDispError     (directInd); 

doMonteCarlo = 'off';

if strcmp(doMonteCarlo,'on')
    
    %% fit through all data
    % make DATA, ERROR and errorModel arrays
    
    DATA            = [allDisp',allPower'];
    DATA(upperBoundInd) = DATA(upperBoundInd)/2;
    
    numData         = length(DATA);
    
    %%%%%%%%% temporary: making error arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dErr                 = allDispError';
    dErr(upperBoundInd)  = allDisp(upperBoundInd)/2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pErr      = log(allPowerError(:,2)) - log(allPower');

    ERROR           = [dErr, pErr];
    
    dispErrorModel  = repmat({'gaussian'},numData,1);
    upperBoundInd   = strcmp(allConstraint,'Upper Bound');
    dispErrorModel(upperBoundInd) ...
                    ={'equal'};
    
    powerErrorModel = repmat({'lognormal'},numData,1);

    errorModel      = [dispErrorModel,powerErrorModel];
    [allDataFit,allDataFitError] = runmontecalrofit(100,DATA,ERROR , ...
                                    'errorModel',      errorModel   , ...
                                    'histogram',       'off'        );
    
    %%  fit through direct data
    % make DATA, ERROR and errorModel arrays
    
    DATA            = [directDisp',directPower'];
    
    numData         = length(DATA);
    
    %%%%%%%%% temporary: making error arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dErr                 = directDispError';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pErr      = log(directPowerErr(:,2)) - log(directPower');

    ERROR           = [dErr, pErr];
    
    dispErrorModel  = repmat({'gaussian'},numData,1);    
    powerErrorModel = repmat({'lognormal'},numData,1);

    errorModel      = [dispErrorModel,powerErrorModel];
    [directDataFit,directDataFitError] = runmontecalrofit(1000,DATA,ERROR   , ...
                                    'errorModel',      errorModel           , ...
                                    'histogram',       'off'                );
                                
    hold on

    monteCarloAllDataFit    = plot(minmaxD, ...
                               10.^(allDataFit(1)*log10(minmaxD)+allDataFit(2))); 

                           
    monteCarloDirectDataFit = plot(minmaxD, ...
                               10.^(directDataFit(1)*log10(minmaxD)+directDataFit(2)));
                            
                                      
end

if strcmp(S.ptest, 'on')
    p = calcSignificance(DATA,directDataFit);
    fprintf('We report a %f condidence of thefit through the directly fit data \n', p)
end

if strcmp(S.histogram, 'on')
    subplot(1,2,1)
end


upperBoundData          = errorbar(upperBoundDisp,  upperBoundPower     , ...
                                  -upperBoundPowerErr(:,1) + upperBoundPower' , ...
                                   upperBoundPowerErr(:,2) - upperBoundPower' , ...
                                   upperBoundDispError                  , ...
                                   upperBoundDispError                  );
                               
directData              = errorbar(directDisp,      directPower         , ...
                                  -directPowerErr(:,1) + directPower'   , ...
                                   directPowerErr(:,2) - directPower'   , ...
                                   directDispError                      , ...
                                   directDispError                      );   

set(upperBoundData, 'LineStyle',        'none'      , ...
                    'Color',            [.7 .7 .7]  , ...
                    'Marker',           's'         , ...
                    'MarkerSize',       5           , ...
                    'MarkerEdgeColor',  [.4 .4 .4]  , ...
                    'MarkerFaceColor',  [1 1 1]     );

set(directData,     'LineStyle',        'none'      , ...
                    'Color',            [.5 .5 .5]  , ...
                    'Marker',           'o'         , ...
                    'MarkerSize',       5           , ...
                    'MarkerEdgeColor',  [0 0 0]     , ...
                    'MarkerFaceColor',  [0 0 0]     );
                
% Zero displacement Data ploted as lines
for iLine = 1:length(zeroDispPower)
    zeroDisplacement = plot(minmaxD, zeroDispPower(iLine)*[1,1]);
    set(zeroDisplacement,   'Color',            [0.5 0.5 0.5], ...
                            'LineWidth',        2            , ...
                            'LineStyle',        '--'         );
end
   
    % pass a best fit line through the entire dataset
    % d = fit(goodData,goodPower,'power1');
    % plot(minmaxD,(d.a*minmaxD.^d.b));
    d = polyfit(log10(allDisp),log10(allPower),1);


allDataFitLine  = plot(minmaxD,10.^(d(1)*log10(minmaxD)+d(2)));

set(allDataFitLine, 'Color',            [1 0 0]      , ...
                    'LineWidth',        2            , ...
                    'LineStyle',        '--'         );         
         
% pass a best fit through well constrained data


% pass a best fit line through the entire dataset
%  d = fit(directDisp ,directPower,'power1','Weigths',(2./diff(directPowerErr)).^2);
%  directDataFit = plot(minmaxD,d(minmaxD), 'Linewidth',2);
 

d = polyfit(log10(directDisp),log10(directPower),1);
directDataFitLine = plot(minmaxD,10.^(d(1)*log10(minmaxD)+d(2)));

set(directDataFitLine,  'Color',            'black'     , ...
                        'Linewidth',        3           );

% graphical considerations

hXLabel = xlabel('Displacement (m)');
hYLabel = ylabel(parameter);
hTitle  = title([parameter, ' as a function of displacement at ',num2str(scale),'m using FFT',' - ' date]);

% hLegend = legend([upperBoundData   , directData, zeroDisplacement   , ...
%                    allDataFitLine, directDataFitLine]               , ...
%                    'Upper bound displacement constraint'            , ...
%                    'Direct displacement constraint'                 , ...
%                    'Zero displacement bound'                        , ...
%                    sprintf('Fit through all data: \\it{P(d) = %0.2g d^{%0.2g \\pm %0.1g}}', ...
%                             10^allDataFit(2), allDataFit(1), allDataFitError(2)),  ...
%                    sprintf('Fit through data with direct displacement constraint: \\it{P(d) = %0.2g d^{%0.2g \\pm %0.1g}}', ...
%                              10^directDataFit(2), directDataFit(1), directDataFitError(2)),  ...
%                    'location','SouthWest'                           );
               


               
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
% set([hLegend, gca]             , ...
%     'FontSize'   , 8           );
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 10          );
set( hTitle                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );       

set(gca,            'XScale',           'log'       , ...
                    'YScale',           'log'       );

set(gca, ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         );
               
% tags to the points for analysis
if strcmp(S.text,'on')
    noLoc       = strcmp(wellProcessed, 'no');
    colorSpec   = zeros(numFiles,3);
    
    for iFile = 1:numFiles
        if noLoc(iFile) == 1
            colorSpec(iFile,:) = [1 0 0];
        end
    end
    offset  = 1.1;
    t       = text(displacement*offset,Power,fileNameArray,...
                  'interpreter', 'none');
    for iFile = 1:numFiles
        t(iFile).Color = colorSpec(iFile,:);
    end
end

hold off

if strcmp(S.histogram, 'on') 
    subplot(1,2,2)
    [counts, bins] = hist(log10(Power));
    barh(bins,counts)
end
end

%% plot Hurst Exponent as a function of Prefactor
    function [] = hurstVsPrefactor(S,files, fileIndex, subsetLoc,numFiles)
    % plots the hurst exponent as a function of theprefacto with the
    % specific intention to extract correlation between these two...
    
orientation = S.orientation;

% added functonality commented out for the moemen (not completed)
% parameter   = S.parameter;     
% scale       = S.scale;   
% displacement= S.displacement;
% constraint  = S.constraint;
% wellProcessed=S.wellProcessed;

hurst       = zeros(1,numFiles);
prefactor   = zeros(1,numFiles);

fileNameArray           = cell(1,numFiles);

for iFile = 1:numFiles
    % get the file and load it
    fileName            = files(fileIndex(subsetLoc(iFile))).name;
    load(fileName,'parameters')
    
     fileNameArray{1,iFile} = parameters.fileName;
    
    % get the roughness data
    spectrumType        = 'FFT';
    desiredData         = parameters.(orientation).(spectrumType);
    fx                  = desiredData{1};
    Px                  = desiredData{1,2};
    Px                  = Px(1:length(fx));
    PxNanInd            = isnan(Px);
    fx                  = fx(~PxNanInd);
    Px                  = Px(~PxNanInd);
    
    % get the fit
    fitObj              = makebestfit(fx,Px,'FitMethod','section','SectionVal',S.fractalSection);
    hurst(iFile)        = (fitObj.p1+1)/-2;
    prefactor(iFile)    = fitObj.p2;
    
    %plot(fitObj,fx,Px)
    %hold on
    
end

% a lot more functionality could be added here...

% make the plot

figure
scatter(prefactor,hurst)
title('Correlation between Hurst exponent and Pre-factor')
xlabel('Prefactor')
ylabel('Hurst Exponent')
hold on

% tags to the points for analysis
offset = 1.5;
text(prefactor * offset, hurst,fileNameArray,...
     'interpreter', 'none');
     
hold off
ax = gca;
set(ax,'XScale', 'log', 'YScale', 'log')

% show statistics of hurst and prefactor

dispStats('hurst',hurst)
dispStats('prefactor',prefactor)

% where

        function dispStats(name, values)
            meanVal       = mean(values);
            rangeVal      = minmax(values);
            formatSpec     = 'The range in values for %s is %f to %f, with a mean value of %f';
            str             = sprintf(formatSpec, name, rangeVal(1), rangeVal(2), meanVal);
            disp(str)
        end
  
figure
histogram(hurst)
    
    end

%% plot a specific, or many, grids (surfaces)
    function [] = plotgrid(inputs)
        
        desiredFileNameArray = inputs{1};
        numGrid  = length(desiredFileNameArray);
        subplotCount = 0;
        
        for iGrid = 1:numGrid
            
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

%% plot all the frequency spectra in a given directory
    function [] = plotspectra(desiredPlot,S,files, fileIndex, subsetLoc,numFiles)
        
   
    % extract all the info from the strucutre now that it is set
    orientation         = S.orientation;          
    displacement        = S.displacement;
    constraint          = S.constraint;
    makeLegend          = S.makeLegend;

    % create plots:
    
    % create the displacement "domain" - for line plotting
    if any(displacement ~= 0)
        minD        = min(displacement(displacement~=0));
        maxD        = max(displacement);
        logMinD = log10(minD);
        logMaxD = log10(maxD);
        logDelD = logMaxD-logMinD;
    end
    
    legendArray = cell(1,numFiles);
    hold on
    
    % loop over files
    for iFile = subsetLoc
        fileName            = files(fileIndex(iFile)).name;
        load(fileName,'parameters')
        
        if strcmp(desiredPlot,'best fits')
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('hi')
            spectrumType        = 'FFT';
            desiredData         = parameters.(orientation).(spectrumType);
            fx                  = desiredData{1,1};
            Px                  = desiredData{1,2};
            
            numFxIn             = length(fx);
            
            Px                  = Px(1:numFxIn);
            
            PxNanInd            = isnan(Px);
            fx                  = fx(~PxNanInd);
            Px                  = Px(~PxNanInd);
            
            try
                errUp               = desiredData{1,3}';
                errDown             = desiredData{1,4}';
                errorArray          = [errUp,errDown];
                errorArray          = errorArray(1:numFxIn,:);
                errorArray          = errorArray(~PxNanInd,:);
            catch
                errorArray          = 'off';
            end

            fitObj              = makebestfit(fx,Px,'FitMethod',     'section'   , ...
                'SectionVal',    0.03        , ...
                'error',         errorArray  );
            plot(fitObj,fx,Px);       
            if ~strcmp(errorArray,'off')
                shadedErrorBar(fx,Px,errorArray','k',1);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             
%             desiredStruct   = getfield(getfield(parameters, orientation),'FFT');
%             desiredCell     = desiredStruct;
%             
%             wavelength      = 1./desiredCell{1};
%             power           = desiredCell{2};
%             
%             power           = power(1:length(wavelength));
%             nanInd          = isnan(power);
%             power           = power(~nanInd);
%             wavelength      = wavelength(~nanInd);
%             
%             d               = polyfit(log10(wavelength'),log10(power),1);
%             x               = [min(wavelength),max(wavelength)];
%             y               = 10.^(d(1)*log10(x)+d(2));
            
        else
            desiredStruct   = parameters.(orientation).(desiredPlot);
            
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
        
        scanName = parameters.fileName;
        legendArray{1,iFile} = scanName; 
            
        % set color as a function of displacement
        if displacement(iFile) == 0
            colorForDisp = [1 0 0];
        else
            colorInd = 0.7*(log10(displacement(iFile))-logMinD)/((logDelD));
            colorForDisp = ones(1,3)-colorInd;
        end
        
        % set unknow displacements to green
        if sum(isnan(colorForDisp)) ~= 0
            colorForDisp = [0, 1, 0];
        end
        
        p = plot(x,y,'-');
        p.Color = colorForDisp;
        formatSpec = 'D = %d - file: %s';
        legendArray{1,iFile} = sprintf(formatSpec, displacement (iFile), ...
            scanName);
        xlabel('Scale (m)')
        
        
        
        % to help distinguish points whith various constraints on displacement
        if strcmp(constraint{iFile}, 'Direct')
            p.LineWidth = 3;
        elseif strcmp(constraint{iFile},'Upper Bound')
            p.LineWidth = 0.2;
        end
        
        
        
        ylabel(desiredPlot)
        title([desiredPlot,' as a function of scale - ', orientation, date])
        
        if strcmp(desiredPlot,'topostd') || strcmp(desiredPlot,'FFT')|| ...
                strcmp(desiredPlot,'PLOMB')
            set(gca,'XScale','log','YScale','log')
            if strcmp(desiredPlot,'FFT')
                ylabel('Power (m^3)')
                xlabel('Wave length, \lambda (m)')                
            else
                xlabel('Scale (log(m))')
            end
        end
    end
    
    
    if strcmp(makeLegend,'yes')
        legend(legendArray,'interpreter', 'none')
    end
    
    end

%% Sub functions
    function [fileNameArray, parameterValueArray, parameterErrorArray] = get_suface_parameter(parameterName, S, files, fileIndex, subsetLoc,numFiles, showOutputYN)
    % function that get a specfic parameter (string: parameterName - one of
    % 'Hurst','prefactor', 'Power','RMS_emily',
    % 'RMS_internet_stylez_with_a_Z'), for a given fileName. Th user also
    % specifies the orientation and type of spectrum (likely 'FFT' that
    % will queried). Note that parameters 'Power','RMS_emily',
    % 'RMS_internet_stylez_with_a_Z' require user to specify a scale of
    % observation
    
    % usage (14-sep-2017):
    % [fileNameArray, Power, powerError] =  get_suface_parameter(fileName,'parallel','FFT','Power'
    
    % author: Kelian Dascher-Cousineau
    % McGill Universityreporting to James Kirkpatric
    
    
    orientation         = S.orientation;
    scale               = S.scale;
      
    % init arrrays
    parameterValueArray         = zeros(1,numFiles);
    parameterErrorArray         = zeros(numFiles,2);
    fileNameArray               = cell(1,numFiles);
    
    for iFile = 1:numFiles
        
        fileName            = files(fileIndex(subsetLoc(iFile))).name;
        load(fileName,'parameters')
        
        fileNameArray{iFile}= parameters.fileName;
        
        % get the roughness data
        spectrumType = 'FFT';
        desiredData         = parameters.(orientation).(spectrumType);
        
        fx                  = desiredData{1,1};
        Px                  = desiredData{1,2};
        
        numFx               = length(fx);
        Px                  = Px(1:numFx);
        
        PxNanInd            = isnan(Px);
        fx                  = fx(~PxNanInd);
        Px                  = Px(~PxNanInd);
        
        % create error bounds or continue with default
        if     strcmp(S.errorBound,'on')
            errUp               = desiredData{1,3};
            errDown             = desiredData{1,4};
            errorBounds         = [errUp, errDown];
            errorBounds         = errorBounds(1:numFx,:);
            errorBounds         = errorBounds(~PxNanInd,:);
        elseif strcmp(S.errorBound,'off')
            errorBounds         = 'off';
        else
            error('error bound must indicate "off" or "on"')
        end
        
        % fit the data
        fitObj                  = makebestfit(fx,Px,'FitMethod',    'section'       , ...
            'SectionVal',   S.fractalSection, ...
            'error',        errorBounds);
        coefConfInt             = confint(fitObj,0.68);
        coefConfInt(:,2)        = 10.^coefConfInt(:,2);
        
        % this is a bit awkward but fuckit im tired (we basically continue
        % with the same code but instead of using the power we just use the
        % parameter specified by the input exponent (the slope of the best
        % fit line through the data).
        
        if      strcmp(parameterName,   'Hurst')
            parameterValueArray(iFile)       	= -(fitObj.p1+1)/2;
            parameterErrorArray(iFile,:)      	= -(coefConfInt(:,1)+1)/2;
            
        elseif  strcmp(parameterName,   'prefactor')
            parameterValueArray(iFile)       	= 10^fitObj.p2;
            parameterErrorArray(iFile,:)      	= coefConfInt(:,2);
            
        elseif  strcmp(parameterName,   'Power')
            parameterValueArray(iFile)        	= 10.^fitObj(log10(1/scale));
            parameterErrorArray(iFile,:)    	= 10.^predint(fitObj,log10(1/scale),0.68);
            
        elseif  strcmp(parameterName,   'RMS_emily')
            % Handle the conversion to RMS as done in Brodsky et al., 2011
            % P(lambda) = C*lambda^BETA
            BETA = -fitObj.p1;
            RMS                     = sqrt(10^fitObj.p2/(BETA - 1))*scale.^((BETA-1)/2);
            parameterValueArray(iFile)          = RMS; %... for the purpose of plotting...
            parameterErrorArray(iFile,:)     	= (coefConfInt(:,2)./(BETAconfInt-1)).^0.5 ...
                .*scale .*(BETAconfInt-1)/2;
            
        elseif strcmp(parameterName,    'RMS_internet_stylez_with_a_Z')
            up2Scale                = fx>1/scale;
            parameterValueArray(iFile)       	= sqrt(trapz(fx(up2Scale),Px(up2Scale)));
            parameterErrorArray(iFile,:)     	= [parameterValueArray(iFile)-10^-20,parameterValueArray(iFile)+10^-20];
        end % if/elseif
    end % for 
    
    if strcmp(showOutputYN, 'yes')
        
        disp([fileNameArray',num2cell(parameterValueArray')])
    elseif ~strcmp(showOutputYN,'no')
        error('Error in function get_suface_parameter input showOutputYN must be ''yes'' or ''no''')
    end
    
    end % function
  
%% small function to query the files in the desired directory
    function [files, numFiles, fileIndex] = getfiles()
    
        directory_name = uigetdir;
        files = dir(directory_name);
        addpath(directory_name)
        
        save('files', 'files')
        fileIndex = find(~[files.isdir]                         & ...
                         ~ strcmp({files.name},'._.DS_Store')   & ...
                         ~ strcmp({files.name},'.DS_Store'));
        numFiles  = length(fileIndex);
    end
    
%% input parsing functions (if you want to add input options here is the place)
    function [S, SUBSETLOC, NUMFILES] = parseInput(inputs, files, numFiles, fileIndex)

% identify data

% possible inputs (not recomended, use addscaninfo instead beforehand)
scanInputInfo = {'displacement'         , ...
                 'constraint'           , ...
                 'displacementError'    , ...
                 'wellProcessed'        , ...
                 'magnification'        , ...
                 'occular'              };
numScanInfo   = length(scanInputInfo);
             
userInput     = {'orientation'          ,...    % slip parallel or perprendicular
                 'parameter'            ,...    % what parameter to plot as a function of displacement
                 'scale'                ,...    % length scale of interpolation
                 'fractalSection'       ,...    % select a section of the spectra to measure slope
                 'subset'               ,...    % select a subset of scans 
                 'text'                 ,...    % make tags next to points
                 'histogram'            ,...    % place histogram of roughness measurements next to plot
                 'errorBound'           ,...    % make error bounds on plot
                 'makeLegend'           ,...    % include legend for spectra plot
                 'bootstrap'            ,...    % run boostrp of the fit
                 'ptest'                };      % evaluate comfidence on fit         
numUserInput  = length(userInput);

numInputs = length(inputs);

if numInputs ~=0
    for iInput = 1:2:length(inputs)
        if ~any([strcmp(inputs(1,iInput),scanInputInfo), ...
                strcmp(inputs(1,iInput),userInput)])
            message = ['input number ',(inputs(1,iInput)),' not allowed'];
            error(message)
        end
    end
end

% default input values:
defaultInput                = [];

defaultInput.orientation    = 'parallel';
defaultInput.parameter      = 'Power';
defaultInput.wellProcssed   = repmat({'yes'},1,numFiles);
defaultInput.constraint     = repmat({'Direct'},1,numFiles);
defaultInput.scale          = 0.01;
defaultInput.fractalSection = 0.03;
defaultInput.text           = 'off';
defaultInput.histogram      = 'off';
defaultInput.errorBound     = 'off';
defaultInput.makeLegend     = 'no';
defaultInput.boostrap       = 'off';
defaultInput.displacementError = ones(1,numFiles); % no great way to do this 
defaultInput.ptest          = 'off';

% query info in the parameter structure of the files

% init arrays
S = defaultInput;
S.constraint                = cell(1,numFiles);
S.displacement              = zeros(1,numFiles);
S.wellProcessed             = cell(1,numFiles);
S.magnification             = zeros(1,numFiles);
S.occular                   = zeros(1,numFiles);

% find information in the inputed directory for each file if it exists
for iFile = 1:numFiles  
    
    fileName            = files(fileIndex(iFile)).name;
    load(fileName,'parameters')
    
    % query fields if they exist
    for iInput = 1:numScanInfo
        f = char(scanInputInfo(1,iInput));
        if isfield(parameters,f)
            if isa(parameters.(f), 'char')
                S.(f){1,iFile}   = ...
                parameters.(f);   
            elseif isa(parameters.(f), 'double')
                S.(f)(1,iFile)   = ...
                parameters.(f); 
            elseif isa(parameters.(f), 'string')
                f = char(f);
                S.(f){1,iFile}   = ...
                parameters.(f); 
            else
                error('unrecognized variable type, must be cell or double')
            end
        end
    end
end

% user specified analysis data (if done so by user)
for iInput = 1:numUserInput
    S       = setVal(S,userInput(1,iInput),inputs);
end

% user specified scan info (if done so by user)
for iInput = 1:numScanInfo
    S       = setVal(S,scanInputInfo(1,iInput),inputs);
end

% select a subset of the files based on a parameter:
if isfield(S,'subset')
    subsetInd   = S.(S.subset{1,1}) == S.subset{1,2};
else 
    subsetInd   = ones(1,numFiles); 
end
SUBSETLOC   = find(subsetInd);
NUMFILES    = sum(subsetInd);

for iScanInfo = 1:numScanInfo          
    S.(scanInputInfo{1,iScanInfo}) = S.(scanInputInfo{1,iScanInfo})(1,SUBSETLOC);
end
   
end


    
