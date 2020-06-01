clear all;

oldpath=path;
path(oldpath,'IMEXfunctions/')
oldpath=path;

global wavenumber mlev which resolution D;

wavenumber=180;

%which IMEX scheme
which=10101;

res=1111; %in NE, 3 km for EC=40e6 m
resolution=res;

%DTset=logspace(log10(0.5),log10(10),10); %works for res=1500
DTset=logspace(log10(1),log10(10),100); %works for res=1500
%DTset=[7];
tt=length(DTset);

%NLEVset=[300];
NLEVset=round(logspace(log10(20),log10(100),50))
NLEVset=unique(NLEVset)

nnl=length(NLEVset);

tol=1e-12;

mmarray=ones(tt,nnl);

for ll=1:nnl
   mlev=NLEVset(ll);

   global_space;
   global_space_imex;
   
   for ii=1:tt
        DT=DTset(ii);
        
        %get Thuburn's matrices
        [nonstiffA,stiffB] = thuburn_system();
        %nonstiffA(:,:)=0;
        chooseIMEX(which);
        
        imexA = IMEXRKstabmat(nonstiffA,stiffB,numvar,DT,A,Ahat,b,bhat,r);
        
        %S1=sparse(imexA);
        kk=eig(imexA);
        %kk=eigs(S1);
        DT
        mlev
        maxeig=max(abs(kk))
        maxeig-1
%         if(maxeig > 1+tol)
%           fillval=2;
%         else
%           fillval=1;
%         end
        mmarray(ii,ll)=maxeig;
   end
end

figure(400000+which);

dzset=D./NLEVset;



tticks=[1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2];
la=1+tticks;

%sort mmarray into bins
mm2=mmarray;
mm2(:,:)=max(max(mmarray));
for i=length(la):-1:1
    la(i)
    nn=find(mmarray<la(i));
    size(nn)
    mm2(nn)=la(i);
end

%now change everything , shift by 1
mm2=mm2-1;
la=la-1;
cp=colormap(parula);
colormap([1 1 1; ...
    0.9769    0.9839    0.0805;...
    0.9628    0.9373    0.1265;...
    0.9786    0.8386    0.1766;...
    0.3406    0.8008    0.4789;...
    0.2440    0.4358    0.9988]);


pcolor(dzset,DTset,mm2);
%shading faceted;
shading flat;

set(gca,'XScale','log');set(gca,'YScale','log')
cb=colorbar
set(gca,'ColorScale','log')
caxis([min(la),max(la)]);
set(cb,'ticks',tticks)

title('$\max(|\lambda|)-1$','Interpreter','latex')

ylabel('$\Delta t$ (sec)','Interpreter','latex')
%xlabel('$n_{\rm {lev}}$','Interpreter','latex')
xlabel('$\Delta z$ (m)','Interpreter','latex')

%set(gca,'XTick',[20 40 57 80 100]);
set(gca,'XTick',[100, 125, 175, 250, 500])
set(gca,'YTick',[1,2,4,8,10,16,20]);


grid on
set(gca,'layer','top')
set(gca, 'XColor', 'black')
set(gca, 'YColor', 'black')
ax=gca;
ax.GridAlpha=1;

set(gca,'FontSize',18)


aa=1


    