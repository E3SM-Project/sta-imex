%stability plot
%Lock's system but Thuburn's plot with dx,dt axes

%clear all;

oldpath=path;
path(oldpath,'IMEXfunctions/')
oldpath=path;
path(oldpath,'plotting/')

%tolerance to define stable eigenvalue
tol=1e-12;
%only dt*kz, dt*kx factors will be used
DT=1;
cr2=1.0045e+05; %copied from thuburn's file
imm=complex(0,1);

%which IMEX scheme
%which=10030;
figure(200000+which);

chooseIMEX(which);

%for loglog spacing
%not really courant number, this is k_x*dt and k_z*dt,
%csound is used in plotting
coux=logspace(-4,.1,nnx);
couz=logspace(-4,2,nnz);


tt=length(couz);
ll=length(coux)

mmarray=ones(tt,ll);
csound = sqrt(cr2);
global_space_imex;

for xind=1:ll
    for zind=1:tt
        
        nonstiffA=-imm*coux(xind)*[0 0 1; 0 0 0; cr2 0 0];
        stiffB   =-imm*couz(zind)*[0 0 0; 0 0 1; 0 cr2 0];
        imexA = IMEXRKstabmat(nonstiffA,stiffB,3,DT,A,Ahat,b,bhat,r);
        
        kk=eig(imexA);
        kk2=sort(abs(kk));
        maxeig=kk2(end)
        
        mmarray(zind,xind) = maxeig;
    end
end

%replace very high values
hi=1.1;
aa=find(mmarray>hi);
mmarray(aa)=hi;

plot_mmarray_pcolor2_courant

set(gca,'FontSize',18)

%recover Ahat last row for title
lastrow=Ahat(end, 1:end);
%parse last row into title
ti=strcat("Courant axes, 2D acoustic, scheme ", int2str(which),", lastrow=[");
for ii=1:length(lastrow)
    ti=strcat(ti,num2str(lastrow(ii)),",");
end
ti=strcat(ti,"]")
%title(ti)
aaa=1