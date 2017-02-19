function [] = gridoverview(zGrid,parameters)
% gives a quick overview of a scan loaded into the directory

numNan = sum(sum(isnan(zGrid)))/numel(zGrid)*100;


figure
subplot(2,2,1)
imagesc(zGrid)
% surf(zGrid*1000000,'EdgeColor','none')
% camlight headlight
% lighting phong
axis equal
colorbar
title('Height field')

subplot(2,2,2)
histogram(zGrid)
titleName = sprintf('Height didstribution, number of nan: %d percent', numNan);
title(titleName)

zpp = del2(zGrid);

subplot(2,2,3)
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




fitObj              = makebestfit(fx,Px,'FitMethod',     'section'   , ...
                                        'SectionVal',    0.03        , ...
                                        'error',         errorArray  );
plot(fitObj,fx,Px)
hold on
shadedErrorBar(fx,Px,errorArray','k',1)
% errorbar(fx, Px,errorArray(:,2),errorArray(:,1))
ax = gca;
set(ax,'XScale', 'log', 'YScale', 'log')
title(['Power Spectrum, Beta =', num2str(fitObj.b)])
xlabel('frequency (m^{-1})')
ylabel('Power (m^3)')

subplot(2,2,4)
histogram(log10(abs(zpp)))
title('Curvature didstribution (log abs)')
end

