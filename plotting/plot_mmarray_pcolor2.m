
plot_prepare1;

courantplot=0;

ressetdx=EC/12./resset/1000;
%making this 2dx wave horizontal axis
%pcolor(ressetdx,DTset,mm2);
pcolor(2*ressetdx,DTset,mm2);
%shading faceted;
shading flat;

plot_prepare2;

ylabel('$\Delta t$ (sec)','Interpreter','latex')
%xlabel('$\Delta x$ (km)','Interpreter','latex')
%xlabel('$2\Delta x = 2\pi/k_x$ (km)','Interpreter','latex')
xlabel('$2\pi/k_x$ (km)','Interpreter','latex')


%set(gca,'XTick',[1,2,4,10,25,50,100])
set(gca,'XTick',2*[1,2,4,10,25,50,100])
set(gca,'YTick',[1,2,4,10,20,50,100,200,300])


grid on
set(gca,'layer','top')
set(gca, 'XColor', 'black')
set(gca, 'YColor', 'black')
ax=gca;
ax.GridAlpha=1;
