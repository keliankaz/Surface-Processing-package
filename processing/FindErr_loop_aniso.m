function [errUp, errDown, Px]=FindErr_loop_aniso(PowerStructx)
% output std for powerspectrum (and mean in case it is lost);

% PSxArray=struct2cell(PowerStructx); % MAKE SURE THIS AGREES WITH THE DIRECTION BEING ANALYZED
% PxiArray_compressed=squeeze(PSxArray(3,:,:));
% clear PSxArray
%  %len = cellfun('prodofsize', PxiArray_compressed);
%   Pxi = zeros(floor(length(PxiArray_compressed)/2), floor(length(PxiArray_compressed{1})));
%   index = 0;
%   maxn=numel(PxiArray_compressed);
%   for n = 1:maxn
%       if (length(isfinite(PxiArray_compressed{n}))>2)
%         index = index + 1;
%         Pxi(index,:) = PxiArray_compressed{n}(:); % more general: c{i}(:)
%       end
%   end
%   % chop off the ends
%   clear PxiArray_compressed;
%   Pxi2=Pxi(1:index,:);
%   clear Pxi;
% %Pxi=cell2mat(PxiArray_compressed')';
% %Errx=10.^nanstd(log10(Pxi));
% Errx=nanstd(Pxi2);
% %Px=10.^nanmean(log10(Pxi2));
% Px=nanmean(Pxi2);
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PSxArray            =struct2cell(PowerStructx); % MAKE SURE THIS AGREES WITH THE DIRECTION BEING ANALYZED
PxiArray_compressed =squeeze(PSxArray(3,:,:));
 
  maxn = numel(PxiArray_compressed);
  Px = zeros(1,maxn);
  errUp = zeros(1,maxn);
  errDown = zeros(1,maxn);
  
  for m = 1:floor(length(PxiArray_compressed{1}))
      [Px(m),errUp(m), errDown(m)] = columnbycolumn(PxiArray_compressed,maxn,m);
  end
end

function [Px, errUp, errDown] = columnbycolumn(PxiArray_compressed,maxn,m)

index = 0;
Pxi = nan(1,maxn);

for n = 1:maxn
    if (length(isfinite(PxiArray_compressed{n}))>2)
        index = index + 1;
        Pxi(n) = PxiArray_compressed{n}(m);
    end
end

Pxi = Pxi(1:index);

%Px = nanmean(Pxi);

% lof transform the data:
logPxi  = log10(Pxi);
logPx   = nanmean(logPxi);
logErr  = nanstd(logPxi);

Px      = 10^logPx;

% asymetrical error bounds
errUp   =   10^(logPx+logErr)   - 10^(logPx); 
errDown = - 10^(logPx-logErr)   + 10^(logPx); 


% take the mode of the data
% NTotal = ceil(sqrt(sum(~isnan(Pxi))));
% if NTotal ~= 0
%     [N,edges] = histcounts(log10(Pxi),NTotal);
%     modeInd = find(N == max(N));
%     modeInd = modeInd(1);
%     Px = 10^((edges(modeInd)+edges(modeInd+1))/2);    
% else
%     Px = nan;
% end

% Errx = nanstd(Pxi);

end