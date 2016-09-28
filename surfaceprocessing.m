function [ ] = surfaceprocessing(unit,toDo,varargin)

% Loops over input xyz files and makes an frequency spectrum.

% Input should be: the unit in which data is imported (i.e. 'mm', 'cm' or 'm') - it is then
% converted into meters

% User will be prompted to select the folder in which the scan files. There
% should be NO other file in the folder. Files should have .xyz format.

% point spacing is determined automatically based on the point density

% output: files will be saved 

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

if nargin >=1
    instrument = varargin{1};
else
    instrument = 'default';
end

if nargin >=2
    numberOfScales = varargin{2};
else
    numberOfScales = 'default'; % default is ten scales
end 

if nargin >=3
    decimationFactor = varargin{3};
else
    decimationFactor = 'default';
end

%% process data:

parfor iFile = 1:length(fileIndex)
    
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
    parameters.NumberOfSampledScales = NumberOfScales;
    parameters.Date         = date;
    
    % save output
    fileNameSpec        = '%s_processing_output.mat';
    outputFileName      = sprintf(fileNameSpec,fileName);
    
    save([destination_directory,'\',outputFileName], ...
        'surface', 'zGrid', 'parameters')
end 
