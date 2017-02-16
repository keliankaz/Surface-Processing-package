function S = setVal(S, name, input)
% set user specified inputs to a default structure based on "pair-wise
% input"

% S         : is a structure to whom 'name' fields will be assigned if they exist in
%             the 'input'
% name      : cell array of all fields to possibly set
% input     : cell array of pair-wise input values

if length(input) >= 2
    
    specifier   = input(1:2:end);
    value       = input(2:2:end);
    
    for iName = 1:length(name)
        
        fieldName   = name(iName);
        ind         = strcmp(fieldName, specifier); % check every odd input
        fieldName   = char(fieldName);
        
        if any(ind)
            S.(fieldName) = value{ind}; %  assign field val to even input
        end
        
    end
end
