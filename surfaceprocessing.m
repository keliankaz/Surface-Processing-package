function [ ] = surfaceprocessing(unit)

% Loops over input xyz files and makes an frequency spectrum.

% Input should be: the unit in which data is imported (i.e. 'mm', 'cm' or 'm') - it is then
% converted into meters

% User will be prompted to select the folder in which the scan files. There
% should be NO other file in the folder. Files should have .xyz format.

% point spacing is determined automatically based on the point density

% output: files will be saved 

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

parfor iFile = 1:length(fileIndex)
    
    fileName            = files(fileIndex(iFile)).name;
    parfor_process(fileName,unit,destination_directory);

end

disp('Alright we are done here!')
end

function [] = parfor_process(fileName,unit,destination_directory)
% this function enssentially enables the parfor loop to be completey
% parallel - otherwise the program runs into transparency issues.

    [surface, zGrid, pointSpacing] = ...
        surface_preprocessing_2(fileName, ...
        unit);
    
    parameters.parallel = ...
        surface_analysis(zGrid,pointSpacing);
    
    parameters.perpendicular = ...
        surface_analysis(zGrid',pointSpacing);
    
    parameters.pointSpacing = pointSpacing;
    parameters.fileName     = fileName;
    
    % save output
    fileNameSpec        = '%s_processing_output.mat';
    outputFileName      = sprintf(fileNameSpec,fileName);
    
    save([destination_directory,'\',outputFileName], ...
        'surface', 'zGrid', 'parameters')
end 
