  function [fitObj] = makebestfit(F,P,varargin)
  % small function to make the best fit of the power spectrum
  %
  %     % ensure the vectors are of the same size/dimension
  %     if length(F) == length(P)
  
  % possible inputs
  inputList                   = { 'FitMethod',... % to fit the data
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
  
  defaultStruct = setVal(defaultStruct,inputList,userSpec);
  
  
  % working structure -
  S           = defaultStruct;
  numIn       = length(F);
  
  if     ~strcmp(S.error,     'off')
      logsigma    = log10(S.error(:,1))-log10(P);
      weightArray = 1./logsigma.^2;
  else
      weightArray = ones(1,numIn);
  end
  
  
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
      
      
      numEntries  = ceil(numIn*S.SectionVal)-1;
      
      f           = F(10:numEntries);
      p           = P(10:numEntries);
      weightArray = weightArray(10:numEntries);
      
  else
      disp('error in the code making the best fit, input vectors must be the same length')
  end
  
  
  fitObj      = fit(log10(f),log10(p),'poly1','Weights',weightArray);
  %             C           = fitObj.p2;
  %             BETA        = fitObj.p1;
  %             H           = (BETA+1)/-2;

    end
