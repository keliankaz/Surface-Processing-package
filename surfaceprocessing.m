function [ ] = surfaceprocessing(varargin)

% Loops over input xyz files and makes an frequency spectrum.

% Inputs:

% all optional, MUST be in pairs with a (the order of pairs is not important): 

% 'unit', followed by the unit in which data is imported (i.e. 'mircon', 'mm', 'cm' or 'm') - it is then
% converted into meters (default is meters)

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

ind = strcmp(varargin,'unit');
if ind ~= 0
    unit = varargin{find(ind)+1};
else
    unit = 'm'; % meters (default)
end 

ind = strcmp(varargin,'toDo');
if ind ~= 0
    toDo = varargin{find(ind)+1};
else
    toDo = 'all'; % meters (default)
end 

% bypass skip preprocessing 'yes', 'no'
ind = strcmp(varargin,'bypass');
if ind ~= 0
    bypass = varargin{find(ind)+1};
else
    bypass = 'no'; % meters (default)
end 

ind = strcmp(varargin,'instrument');
if ind ~= 0
    instrument = varargin{find(ind)+1};
else
    instrument = 'default';
end

ind = strcmp(varargin,'numberOfScales');
if ind ~= 0
    numberOfScales = varargin{find(ind)+1};
else
    numberOfScales = 'default'; % default is ten scales
end 

ind = strcmp(varargin,'decimationFactor');
if ind ~= 0
    decimationFactor = varargin{find(ind)+1};
else
    decimationFactor = 'default';
end

%% process data:

for iFile = 1:length(fileIndex)
    
    fileName            = files(fileIndex(iFile)).name;
    parfor_process(fileName,unit,toDo,destination_directory,instrument,...
                   numberOfScales, decimationFactor);
    
end

%%
disp('Alright we are done here!')
end

function [] = parfor_process(fileName,unit,toDo,destination_directory, ...
                             instrument,numberOfScales,decimationFactor)
% this function enssentially enables the parfor loop to be completey
% parallel - otherwise the program runs into transparency issues.



    [surface, zGrid, pointSpacing] = ...
        surface_preprocessing_2(fileName,unit,instrument);
        
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
    
    % save output
    fileNameSpec        = '%s_processing_output.mat';
    outputFileName      = sprintf(fileNameSpec,fileName);
    
    save([destination_directory,'\',outputFileName], ...
        'surface', 'zGrid', 'parameters')
end 