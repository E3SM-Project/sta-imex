%Lock's system but Thuburn's plot

clear all;

oldpath=path;
path(oldpath,'IMEXfunctions/')
oldpath=path;
path(oldpath,'plotting/')

%global vars defined here before calling global_space*:
global wavenumber which resolution mlev imm k cr2 D;

%wavenumber, num of vert levels
wavenumber=180;
mlev=72;

%which IMEX scheme
which=10030;

figure(10000+which);

%production :16, :0.5
%mid res :60, 1
%coarse :120, 2
%set of hor resolutions, in NE
%resset=[30:60:3000]; %high res 3000,
%to compute set of dt
%ibdt=[1:1:20];  %try 0.25
%set of dt
%DTset=1/2* (1.4).^ibdt ;


resset=logspace(log10(30),log10(3500),100);
DTset=logspace(log10(0.5),log10(400),100);
ll=length(resset);
tt=length(DTset);

%tolerance to define stable eigenvalue
tol=1e-12;
chooseIMEX(which);
%time and resolution array

%one set of kz
%dD=10^4;
%dzset = [2:mlev].*dD./mlev;
%newlev=length(dzset);

%another set of kz -- based on Lock's kz/kx \in [1e-2, 1e4]
kzset=logspace(-2,4);
newlev=length(kzset);

mmarray=ones(tt,ll);
for resind=1:ll
    resolution=resset(resind);
    global_space;
    %c_sound, D defined in glob space
    csound = sqrt(cr2);
    global_space_imex;
    
    for ii=1:tt
        DT=DTset(ii);
        
        %vert wavenumber array
        kzarray=ones(newlev,1);
        for nlev=1:newlev
            
            %Tlength=dzset(nlev);
            %kz=2*pi/Tlength;
            
            kz=kzset(nlev)*k;
            
            nonstiffA=-imm*k *[0 0 1; 0 0 0; cr2 0 0];
            stiffB   =-imm*kz*[0 0 0; 0 0 1; 0 cr2 0];
            %stiffB(:,:) = 0;
            imexA = IMEXRKstabmat(nonstiffA,stiffB,3,DT,A,Ahat,b,bhat,r);
            
            kk=eig(imexA);
            kk2=sort(abs(kk));
            maxeig=kk2(end);
            
            %[resolution DT maxeig]
            kzarray(nlev)=maxeig;
            
        end
       
        DT 
        max(kzarray)
        
        mmarray(ii,resind) = max(kzarray);
    end
    
    resolution
end

%replace very high values
hi=1.1;
aa=find(mmarray>hi);
mmarray(aa)=hi;

plot_mmarray_pcolor2

set(gca,'FontSize',18)

%recover Ahat last row for title
lastrow=Ahat(end, 1:end);
%parse last row into title
ti=strcat("2D acoustic, scheme ", int2str(which),", lastrow=[");
for ii=1:length(lastrow)
    ti=strcat(ti,num2str(lastrow(ii)),",");
end
ti=strcat(ti,"]")
%title(ti)
aaa=1