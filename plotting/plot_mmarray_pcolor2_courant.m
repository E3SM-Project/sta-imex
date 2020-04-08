
plot_prepare1;

courantplot=1;

%plot with dx
pcolor(sqrt(cr2)*coux,sqrt(cr2)*couz,mm2);
shading flat;

plot_prepare2;

ylabel('$k_z \Delta t$ (sec/m)','Interpreter','latex')
xlabel('$k_x \Delta t$ (sec/m)','Interpreter','latex')

set(gca,'XTick',[0.1,1,2,3,4])
set(gca,'YTick',[0.1,1,2,4])
grid on
set(gca,'layer','top')
set(gca, 'XColor', 'black')
set(gca, 'YColor', 'black')
ax=gca;
ax.GridAlpha=1;