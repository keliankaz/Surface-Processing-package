function [ ] = surfaceprocessing(varargin)

% Loops over input files of spactial data and performs desired analyses
% including FFT and PLOMB spectra; RMS , kurtosis, skewness parameters; ALL
% along  profiles. Preprocessing can be performed to remove defects,
% rotate data and align grids. 

% Inputs:

% all optional, MUST be in pairs with a (the order of pairs is not important): 

% 'unit', followed by the unit in which data is imported (i.e. 'mircon',
% 'mm', 'cm' or 'm') - it is then converted into meters (default is meters)
% - do not need to specify if bypass for the zygo format is activated

% 'toDo', followed by the desired analyses on of: 'FFT','PLOMB',
% 'parameters' or 'all' (default is 'all') - can be a cell array

% 'bypass', followed by 'zygo' 'pre-processing' or 'no' to  be used
% input is already in aligned clean grid form - input files are then
% (default is 'no'). zygo is adapted to the proprietary data format of the
% white light in wong. the 'pre-processing' option requires a .mat
% structure with a field named 'grid' with the topography and a field name
% 'pointSpacing' specifying the point spacing (in meters). In either case
% the topography must be aligned such that the positive x direction is the
% parallel direction

% 'numberOfScales' followed by the desired number of analysed scales. THis
% has a lot of effect on the amount of processing time (default is 10)

% 'decimationFactor' followed by the desired decimation factor (default is
% 1)

% 'instrument' followed by 'white light', 'laser scanner' or 'lidar'
% (default does not set any instrument specific adjustments

% User will be prompted to select the folder in which the scan files. There
% should be NO other file in the folder. Files should have .xyz format.

% point spacing is determined automatically based on the point density

% output: files will be saved 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ex.

% surfaceprocessing('bypass', 'zygo','toDo', 'FFT','decimationFactor',5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% running surface asymmetry function and saving files

tic

disp('Choose the input directory')
directory_name = uigetdir;
disp('good job!')
files = dir(directory_name);
addpath(directory_name)

disp('Now choose the output directory')
fileIndex = find(~[files.isdir] & ~strcmp({files.name},'.DS_Store'));

destination_directory = uigetdir;
disp('Nice!')

%% dealing with user inputs

% possible user inputs:

userInputs      = {'unit',...
                   'toDo',...
                   'bypass',...
                   'instrument',...
                   'numberOfScales',...
                   'decimationFactor'};

% default values;
default                 = [];       % initialize structure
default.unit            = 'm';      % meters
default.toDo            = 'all';    % runs all analyses
default.bypass          = 'no';     % does not by pass pre-processing step
default.instrument      = 'default';% runs without any specific instrument preferences
default.numberOfScales  = 20;       % the number of scals(log space analysed
default.decimationFactor= 1;        % decimation of point cloud

S = setVal(default,userInputs,varargin);

% % unit
% ind = strcmp(varargin,'unit');
% if any(ind)
%     unit = varargin{find(ind)+1};
% else
%     unit = 'm'; % meters (default)
% end 
% 
% % to do list
% ind = strcmp(varargin,'toDo');
% if any(ind)
%     toDo = varargin{find(ind)+1};
% else
%     toDo = 'all'; % meters (default)
% end 
% 
% % bypass skip preprocessing 'yes', 'no'
% ind = strcmp(varargin,'bypass');
% if any(ind)
%     bypass = varargin{find(ind)+1};
% else
%     bypass = 'no'; % (default)
% end 
% 
% % instrument
% ind = strcmp(varargin,'instrument');
% if any(ind)
%     instrument = varargin{find(ind)+1};
% else
%     instrument = 'default';
% end
% 
% % number of scale
% ind = strcmp(varargin,'numberOfScales');
% if any(ind)
%     numberOfScales = varargin{find(ind)+1};
% else
%     numberOfScales = 'default'; % default is ten scales
% end 
% 
% % decimation factor
% ind = strcmp(varargin,'decimationFactor');
% if any(ind)
%     decimationFactor = varargin{find(ind)+1};
% else
%     decimationFactor = 'default';
% end

%% process data:

parfor iFile = 1:length(fileIndex)
    tic 
    fileName            = files(fileIndex(iFile)).name;
    disp(['Now processing: ',fileName])
    parfor_process(fileName,S.unit,S.toDo,destination_directory, S.bypass, ...
                   S.instrument, S.numberOfScales, S.decimationFactor);
    oneFileTime         = toc;
    
    disp(['last file took', num2srt(toc),'seconds'])
    time2finish         = oneFileTime*(length(fileIndex)-iFile);
    disp(['estimated time to finish: ',num2str(time2finish)]);
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
        
    elseif strcmp(bypass,'zygo')
        % import, parse and detrend data
        [zGrid,pointSpacing] = parse_zygo_format('fileName',fileName, ...
                                                 'detrend','yes');
        surface     = zGrid;
        
    elseif strcmp(bypass,'pre-processing')
        structIn    = load(fileName);
        zGrid       = structIn.grid;
        pointSpacing= structIn.pointSpacing;
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