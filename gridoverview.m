function [] = gridoverview(zGrid)
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
imagesc(log10(abs(zpp)))
colorbar
axis equal
title('log curvature field')

subplot(2,2,4)
histogram(log10(abs(zpp)))
title('Curvature didstribution (log abs)')
end