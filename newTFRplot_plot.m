
function []=newTFRplot_plot(T, F, D, S, z_lim, patchgirth,framewidth, lettersize, LIM,...
    da_yticks,xLIM,col_map,plot_title)
[~, daLIM1]=min(abs(S-LIM(1)));daLIM1=daLIM1(1);
[~, daLIM2]=min(abs(S-LIM(2)));daLIM2=daLIM2(1);
for iter=1:length(da_yticks)
[~, dummy]=min(abs(S-da_yticks(iter)));daT(iter)=dummy(1);
end
imagesc(T, F, D, z_lim);
for iter=1:length(patchgirth)
patch([patchgirth{iter}(1) patchgirth{iter}(2) patchgirth{iter}(2) patchgirth{iter}(1)],[daLIM1 daLIM1 daLIM2 daLIM2], 'black');
end
ax = gca;
ax.LineWidth = framewidth;
ax.FontSize = lettersize;
ax.FontWeight = 'bold';
yticks(daT)
yticklabels(cellfun(@num2str, num2cell(da_yticks)','un', 0))
xlim(xLIM)
ylim([daLIM1 daLIM2])
colormap(col_map)
set(gca,'YDir','normal')
title(plot_title)
end







