function [] = gridoverview(zGrid,parameters)
% gives a quick overview of a scan loaded into the directory

numNan = sum(sum(isnan(zGrid)))/numel(zGrid)*100;


figure
subplot(2,2,1)
surf(zGrid*1000000,'EdgeColor','none')
camlight headlight
lighting phong
axis equal
colorbar
title('Height field (100000000x vert exageration)')

subplot(2,2,2)
histogram(zGrid)
titleName = sprintf('Height didstribution, number of nan: %d percent', numNan)
title(titleName)

zpp = del2(zGrid);

subplot(2,2,3)
% imagesc(log10(abs(zpp)))
% colorbar
% axis equal
% title('log curvature field')
spectrumType        = 'FFT';
desiredData         = parameters.parallel.(spectrumType);
fx                  = desiredData{1};
Px                  = desiredData{1,2};
Px                  = Px(1:length(fx));
PxNanInd            = isnan(Px);
fx                  = fx(~PxNanInd);
Px                  = Px(~PxNanInd);
plot(fx,Px,'.-')
hold on
[C,H]               = makebestfit(fx,Px,'FitMethod','section','SectionVal',0.03);
dataFit = C*fx.^(-1-2*H);
plot(fx,dataFit, 'r')
ax = gca;
set(ax,'XScale', 'log', 'YScale', 'log')
title(['Power Spectrum, H =', num2str(H)])
xlabel('frequency (m^{-1})')
ylabel('Power (m^3)')

subplot(2,2,4)
histogram(log10(abs(zpp)))
title('Curvature didstribution (log abs)')
end

