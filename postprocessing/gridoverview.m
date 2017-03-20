function [] = gridoverview(zGrid,parameters)
% gives a quick overview of a scan loaded into the directory

numNan = sum(sum(isnan(zGrid)))/numel(zGrid)*100;
verticalExageration = 0.1;
[numY, numX]        = size(zGrid);
y                   = 1:numY;
x                   = 1:numX;
[X, Y]              = meshgrid(x,y);
xGrid               = X*parameters.pointSpacing;
yGrid               = Y*parameters.pointSpacing;
yLim1 = prctile(zGrid(:),0.1);
yLim2 = prctile(zGrid(:),99.9);

figure

% set the location of the plot that will be displayed in figure window,
% location is normalized to 1

ax1 = axes( 'Position',         [0.1  0.1 0.35 0.8 ]);
ax2 = axes('Position',          [0.35 0.7 0.1  0.15]);
ax3 = axes('Position',          [0.55 0.1 0.35 0.8 ]);

faultSurf       = surf(ax1,xGrid,yGrid,zGrid);
faultSurfLight  = light;
faultSurfCBar = colorbar(ax1);


set(faultSurf,      'EdgeColor',    'none');
set(faultSurfLight, 'Parent',           ax1   );
set(ax1,            'ZLimMode',         'manual'                        ,...
                    'ZLim',             [nanmin(nanmin(zGrid))*4, nanmax(nanmax(zGrid))*4], ...
                    'XlimMode',         'manual'                        ,...
                    'XLim',             [nanmin(nanmin(xGrid))  , nanmax(nanmax(xGrid))  ], ...                   ,...
                    'YLim',             [nanmin(nanmin(yGrid))  , nanmax(nanmax(yGrid))  ]);
set(faultSurfCBar,  'Box',              'off'                           ,...
                    'location',         'southoutside'                  ,...
                    'limits',           [yLim1,yLim2]                   );
faultSurfCBar.Label.String = 'm';

                    
                    
                
camlight(faultSurfLight, 'right')
lighting('phong')


zlabel('height (m)')
title(ax1, 'Height field')

% maxVal = max(max(zGrid));
% minVal = min(min(zGrid));
% 
% zGridScaled = (zGrid-minVal)/maxVal;
% zGridSize = size(zGrid);
% colorMap = ones(zGridSize(1),zGridSize(2),3).*zGridScaled;
% 
% zp = gradient(zGrid);
% maxZp = max(max(zp));
% minZp = min(min(zp));
% zpScaled = (zp-minZp)/maxZp;
% gradMap = ones(zGridSize(1),zGridSize(2),3).*zpScaled;
% 
% colorMap = colorMap.*gradMap;

% figure
% image(colorMap)

% imagesc(ax1, zGrid)
% axis equal
% colorbar
% title('Height field')

[N,edges]   = histcounts(zGrid,500);
centers = mean([edges(1:end-1);edges(2:end)]);
barh(ax2,centers,N)
titleName = sprintf('Height didstribution from mean plane (n = %d)', numNan);
title(titleName)
xlabel('Height (m)')
set(ax2,            'XTickLabel',   []    ,...
                    'YLim',         [yLim1,yLim2]);
                    

% zpp = del2(zGrid);

% imagesc(log10(abs(zpp)))
% colorbar
% axis equal
% title('log curvature field')

spectrumType        = 'FFT';
desiredData         = parameters.parallel.(spectrumType);
fx                  = desiredData{1,1};
Px                  = desiredData{1,2};

numFxIn             = length(fx);

Px                  = Px(1:numFxIn);

PxNanInd            = isnan(Px);
fx                  = fx(~PxNanInd);
Px                  = Px(~PxNanInd);

 try 
    errUp               = desiredData{1,3}';
    errDown             = desiredData{1,4}';
    errorArray          = [errUp,errDown];
    errorArray          = errorArray(1:numFxIn,:);
    errorArray          = errorArray(~PxNanInd,:);
 catch
     errorArray          = 'off';
 end



sectionVal         = 0.2;
fitObj              = makebestfit(fx,Px,'FitMethod',     'section'   , ...
                                        'SectionVal',    sectionVal  , ...
                                        'error',         errorArray  );
fractalFitLine  = plot(fx,10.^(fitObj.p1*log10(fx)+fitObj.p2),'r');
hold on

if ~strcmp(errorArray,'off')
    shadedErrorBar(fx,Px,errorArray','k',1);
end

rejectStart     = ceil(length(fx)*sectionVal);
dataScatterGood = plot(fx(1:rejectStart)    , Px(1:rejectStart)     , 'b.-');
dataScatterBad  = plot(fx(rejectStart:end)  , Px(rejectStart:end)   , 'r.-');

level           = 0.95;
xConfInt        = logspace(log10(min(fx)),log10(max(fx)), 50);
yConfInt        = predint(fitObj, log10(xConfInt), level);
confIntLines    = plot(ax3, xConfInt,10.^(yConfInt),'m--');

legend([dataScatterGood, dataScatterBad, fractalFitLine], ...
        'power spectral density (geometrical mean +/- 1\sigma)' , ...
        'rejected data'                                         , ...          
        sprintf('fit through fractal bands width: p(f) = %0.1d f^{%0.1f}',10^fitObj.p2, fitObj.p1));


ax = gca;
set(ax,'XScale', 'log', 'YScale', 'log')
title('Power Spectrum')
xlabel('frequency (m^{-1})')
ylabel('Power (m^3)')

        
% subplot(1,4,4)
% histogram(log10(abs(zpp)))
% title('Curvature didstribution (log abs)')
end

