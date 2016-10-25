function [ ] = surfaceprocessing(varargin)

% Loops over input xyz files and makes an frequency spectrum.

% Inputs:

% all optional, MUST be in pairs with a (the order of pairs is not important): 

% 'unit', followed by the unit in which data is imported (i.e. 'mircon', 'mm', 'cm' or 'm') - it is then
% converted into meters (default is meters) - do not need to specify if
% bypass for the zygo format is activated

% 'toDo', followed by the desired analyses on of: 'FFT','PLOMB',
% 'parameters' or 'all' (default is 'all') - can be a cell array

% 'bypass', followed by 'yes' or 'no' to  be used input is already in aligned clean grid form - input
% files are then (default is 'no'). At the moment bypass is only available
% for the proprietary white light format which includes point spacing
% information in the file. Bypass activates the file parsing function to do
% so.

% 'numberOfScales' followed by the desired number of analysed scales. THis
% has a lot of effect on the amount of processing time (default is 10)

% 'decimationFactor' followed by the desired decimation factor (default is
% 1)

% 'instrument' followed by 'white light', 'laser scanner' or 'lidar'
% (default does not set any instrument specific adjustments

% User will be prompted to select the folder in which the scan files. There
% should be NO other file in the folder. Files should have .xyz format.

% point spacing is determined automatically based on the point density

% todo cell array specifies desired analyses to do: 

% bypass - 

% output: files will be saved 

%%%%%%%%%%%%%%
%ex.

% surfaceprocessing('unit', 'mm','instrument','white light','bypass','yes')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% running surface asymmetry function and saving files

tic

disp('Choose the input directory')
directory_name = uigetdir;
disp('good job!')
files = dir(directory_name);
addpath(directory_name)

disp('Now choose the output directory')
fileIndex = find(~[files.isdir]);
disp('Nice!')

destination_directory = uigetdir;

%% dealing with user inputs

% unit
ind = strcmp(varargin,'unit');
if any(ind)
    unit = varargin{find(ind)+1};
else
    unit = 'm'; % meters (default)
end 

% to do list
ind = strcmp(varargin,'toDo');
if any(ind)
    toDo = varargin{find(ind)+1};
else
    toDo = 'all'; % meters (default)
end 

% bypass skip preprocessing 'yes', 'no'
ind = strcmp(varargin,'bypass');
if any(ind)
    bypass = varargin{find(ind)+1};
else
    bypass = 'no'; % (default)
end 

% instrument
ind = strcmp(varargin,'instrument');
if any(ind)
    instrument = varargin{find(ind)+1};
else
    instrument = 'default';
end

% number of scale
ind = strcmp(varargin,'numberOfScales');
if any(ind)
    numberOfScales = varargin{find(ind)+1};
else
    numberOfScales = 'default'; % default is ten scales
end 

% decimation factor
ind = strcmp(varargin,'decimationFactor');
if any(ind)
    decimationFactor = varargin{find(ind)+1};
else
    decimationFactor = 'default';
end

%% process data:

for iFile = 1:length(fileIndex)
    
    fileName            = files(fileIndex(iFile)).name;
    parfor_process(fileName,unit,toDo,destination_directory, bypass, ...
                   instrument, numberOfScales, decimationFactor);
    
end

%%
disp('Alright we are done here!')
end

function [] = parfor_process(fileName,unit,toDo,destination_directory, ...
                             bypass,instrument,numberOfScales, ...
                             decimationFactor)
% this function enssentially enables the parfor loop to be completey
% parallel - otherwise the program runs into transparency issues.

    if strcmp(bypass,'no')
        [surface, zGrid, pointSpacing] = ...
            surface_preprocessing_2(fileName,unit,instrument);
        
    elseif strcmp(bypass,'yes')
        % import, parse and detrend data
        [zGrid,pointSpacing] = parse_zygo_format('fileName',fileName, ...
                                                 'detrend','yes');
        surface = zGrid;
    else
        disp('warning bypass must be yes or no')
    end
        
    parameters.parallel = ...
        surface_analysis(zGrid,pointSpacing,numberOfScales, ...
                         decimationFactor,toDo);
    
    parameters.perpendicular = ...
        surface_analysis(zGrid',pointSpacing,numberOfScales, ...
                         decimationFactor,toDo);
    
    parameters.pointSpacing = pointSpacing;
    parameters.fileName     = fileName;
    parameters.Instrument   = instrument;
    parameters.Decimation   = decimationFactor;
    parameters.NumberOfSampledScales = numberOfScales;
    parameters.Date         = date;
    parameters.processingTime = toc;
    
    % save output
    fileNameSpec        = '%s_processing_output.mat';
    outputFileName      = sprintf(fileNameSpec,fileName);
    
    save([destination_directory,'\',outputFileName], ...
        'surface', 'zGrid', 'parameters')
end 