%plot IMEx stability plot using dx,dt for axes, our most used plot
%clear all;

oldpath=path;
path(oldpath,'IMEXfunctions/')
oldpath=path;
path(oldpath,'plotting/')

%global vars defined here before calling global_space*:
global wavenumber mlev which resolution;

%wavenumber based on resolution is computed in global_space
%wavenumber const, num of vert levels
wavenumber=180;
mlev=72;

%which IMEX scheme
%which=10103;

figure(20000000+which);

 
%%%%%%%% PAPER, resset is in NE, homme notations for hor. res.
resset=logspace(log10(30),log10(3500),nnx); %high res 3000,
%coarse
%resset=logspace(log10(30),log10(3500),30); %high res 3000,
ll=length(resset);

%to compute set of dt
%ibdt=[1:1:20];  %try 0.25
%DTset=1/2* (1.4).^ibdt ;

%%%%%%%%%PAPER
DTset=logspace(log10(0.5),log10(400),nnz);
%coarse
%DTset=logspace(log10(0.5),log10(400),30);
tt=length(DTset);

%tolerance to define stable eigenvalue
tol=1e-12;

mmarray=ones(tt,ll);

for resind=1:ll
    resolution=resset(resind);
    global_space;
    global_space_imex;
    for ii=1:tt
        DT=DTset(ii);
        
        %get Thuburn's matrices
        [nonstiffA,stiffB] = thuburn_system();
        
        chooseIMEX(which);
        
        imexA = IMEXRKstabmat(nonstiffA,stiffB,numvar,DT,A,Ahat,b,bhat,r);
        
        kk=eig(imexA);
        maxeig=max(abs(kk))
        
        resolution 
        DT 
        maxeig
        mmarray(ii,resind)=maxeig;
        
    end
end

%replace very high values
hi=1.1;
aa=find(mmarray>hi);
mmarray(aa)=hi;

plot_mmarray_pcolor2
set(gca,'FontSize',18)
 
% % %recover Ahat last row for title
% lastrow=Ahat(end, 1:end);
% % %parse last row into title
% ti=strcat("5-waves, scheme ", int2str(which),", lastrow=[");
% for ii=1:length(lastrow)
%     ti=strcat(ti,num2str(lastrow(ii)),",");
% end
% ti=strcat(ti,"]")
% title(ti)
aaa=1





