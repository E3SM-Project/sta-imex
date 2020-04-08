clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% for these use 20 levels only and
%%%%%%%%%%%%%%%%%%%%%%%%%%%% original k as in dispersion plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%% and in Thuburn's paper, k=2*pi/10^6


oldpath=path;
path(oldpath,'IMEXfunctions/')
oldpath=path;
path(oldpath,'plotting/')

%global vars defined here before calling global_space*:
global wavenumber mlev which resolution;

resolution=1;
wavenumber=360*40/12;

mlev=20;

%which IMEX scheme
which=10096;

DT=50;

global_space;
global_space_imex;
[nonstiffA,stiffB] = thuburn_system();
        
chooseIMEX(which);
        
imexA = IMEXRKstabmat(nonstiffA,stiffB,numvar,DT,A,Ahat,b,bhat,r);
[imexAvec,imexAeig]=eig(imexA);

ti=strcat('for-spacetime-dispersion/method-mlev',int2str(mlev),'-',int2str(which),'-dt',num2str(DT),'.mat');
save(ti,'imexA','imexAvec','imexAeig','DT');




