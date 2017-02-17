  function [C,H] = makebestfit(F,P,varargin)
        % small function to make the best fit of the power spectrum
        %
        %     % ensure the vectors are of the same size/dimension
        %     if length(F) == length(P)
        
        % possible inputs
        inputList                   = {'FitMethod',... % to fit the data
                                       'sectionVal',...% section to be fit (0 to 1)
                                       'avgWindow',... % moving average windoe
                                       'smoothing',... % smooth data ('yes'/'no')
                                       'rolloff',...   % starting point of the section
                                       'error'};       % error bounds on data 
                                   
        numInputList                = length(inputList);
        
        % make default setting ( caution, values inputed into function
        % override these values )
        defaultStruct.FitMethod     = 'section';
        defaultStruct.SectionVal    = 0.2;
        defaultStruct.avgWindow     = 11;
        defaultStruct.smoothing     = 'yes';
        defaultStruct.rolloff       = 10;
        defaultStruct.error         = 'off';    % equal weight to all values
        
        userSpec                    = varargin;
        
        % user specified inputs
        for iInput = 1:numInputList
            defaultStruct = setVal(defaultStruct,inputList(1,iInput),userSpec);
        end
        
        % working structure - 
        S           = defaultStruct;
        
        if      size(F) ~= size(P);          P = P';                     end
        if      strcmp(S.smoothing, 'yes');  P = movmean(P,S.avgWindow); end
        
        if      strcmp(S.FitMethod, 'default')
            % Simplest form: take the fit over the entire calculated spectrum:
            f           = F;
            p           = P;
            
        elseif  strcmp(S.FitMethod, 'section')
            % Uninspired, but better form: take the best fit only on a
            % specifict section of the fit.
            
            % In the case of the scan data, this
            % only takes into account the larger wave lengths where the data
            % should not be affected by instrumental artefacts
            
            numIn       = length(F);
            numEntries  = ceil(numIn*S.SectionVal)-1;
            
            f           = F(10:numEntries);
            p           = P(10:numEntries);

        else
            disp('error in the code making the best fit, input vectors must be the same length')
        end
        
        if     ~strcmp(S.error,     'off')
            weigths = log10(1/(S.error(:,1)-S.error(:,2)).^2)';

            fitObj      = fit(f,p,'power1');
            C           = fitObj.a;
            BETA        = fitObj.b;
            H           = (BETA+1)/-2;

    end
