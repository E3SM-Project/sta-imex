set(gca,'XScale','log');set(gca,'YScale','log')
%cb=colorbar('YTick',la);cb.Ruler.Scale='log';cb.Ruler.MinorTick='on';
cb=colorbar
set(gca,'ColorScale','log')
caxis([min(la),max(la)]);
set(cb,'ticks',tticks)

title('$\max(|\lambda|)-1$','Interpreter','latex')