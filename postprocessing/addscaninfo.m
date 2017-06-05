function addscaninfo(varargin)
% function to set additional information to scans

% input: in pair wise form

% unlimited amount of pairs
% 1- a string specifying the field that will be created/updated
% 2- cell or double array with value of the fields. Array must be the same
% length as the number of files that are updated

% user is promted to navigate to desired directory. Thedirectory can only
% have the matlab workspaces of the scans

% CAUTION: THIS OVERWRITES THE ORIGINAL FILES - BE CAREFUL

[files, numFiles,fileIndex] = getfiles();

% loop over files
for iFile = 1:numFiles
    
    fileName                = files(fileIndex(iFile)).name;
    load(fileName,'parameters')
    fileNameArray{1,iFile}  = getfield(parameters,'fileName');
    
    % loop over input pairs
    for iItem = 1:2:(length(varargin)-1)
        
        fieldName               = varargin{1,iItem};
        parameters              = setval(parameters,fieldName,varargin);

    end
    
    save([directory_name,'\',fileName], 'parameters', '-append')
end

    function updatedS = setval(S, name, input)
        % set user specified inputs to a default structure based on "pair-wise
        % input"
        ind = strcmp(name, input);
        value = input{1,find(ind)+1};
        
        if isa(value,'cell') 
            S.(name) = value{iFile};
        elseif isa(value,'double')
            S.(name) = value(iFile);
        else
            error('Input %f must be double or cell array.',iItem)
        end
        
        updatedS = S;
    end

    function [files, numFiles, fileIndex] = getfiles()
        directory_name = uigetdir;
        files = dir(directory_name);
        addpath(directory_name)
        save('files', 'files')
        
        fileIndex = find(~[files.isdir]& ~ strcmp({files.name},'._.DS_Store') & ~ strcmp({files.name},'.DS_Store'));
        numFiles  = length(fileIndex);
    end

end




