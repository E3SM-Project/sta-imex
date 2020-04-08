clear all;

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

figure(300001);
%plotting
LL=2.0;
MS=10;

% ~ 1 km 3072
resolution=3072; %high res 3000,
DTset= [250:2:400] * 30/resolution/1.5; %ne30 set, scale it
tt=length(DTset);

tol=1e-12;
%which IMEX scheme

global_space;
global_space_imex;


whichset=[10090 10091 10100 10101 10102];
whichlabel={"ttype9","ttype9mod","ttype10","ttype10mod","ttype10as"};
whichcurve=["b-","black-","green-","r-","m-"];

if((length(whichset) ~= length(whichlabel)) || (length(whichset) ~= length(whichcurve)))
    stop
end

ww=length(whichset);
mmarray=ones(tt,ww);


global cr2 imm k;

%another set of kz -- based on Lock's kz/kx \in [1e-2, 1e4]
kzset=logspace(-2,4);
newlev=length(kzset);

for jj=1:ww
    
    which = whichset(jj);
    
    chooseIMEX(which);
    
    for ii=1:tt
        
        DT=DTset(ii);
        %vert wavenumber array
        kzarray=ones(newlev,1);
        for nlev=1:newlev
            kz=kzset(nlev)*k;
            nonstiffA=-imm*k *[0 0 1; 0 0 0; cr2 0 0];
            stiffB   =-imm*kz*[0 0 0; 0 0 1; 0 cr2 0];
            
            imexA = IMEXRKstabmat(nonstiffA,stiffB,3,DT,A,Ahat,b,bhat,r);
            kk=eig(imexA);
            mmax=max(abs(kk));
            kzarray(nlev)=mmax;
        end
        
        mmax=max(kzarray);
        %define all next dt unstable if this one is unstable
        if(mmax > 1+tol)
            mmarray(ii:tt,jj) = 2;
            disp("exiting")
            which
            break
        end
        
    end
    
    %plot(DTset,mmarray(:,jj),whichcurve(jj),'DisplayName',whichlabel(jj)); hold on;
    plot(DTset,mmarray(:,jj),whichcurve(jj),'LineWidth',LL); hold on;
    
    which
end
set(gca,'FontSize',18)
axis([DTset(1) DTset(end) 0.5 2.5])
legend(whichlabel,'Location','northwest','NumColumns',2)
ti=strcat("res ",int2str(resolution),"ne")
title(ti)

% %recover Ahat last row for title
% lastrow=Ahat(end, 1:end);
% %parse last row into title
% ti=strcat("5-waves, scheme ", int2str(which),", lastrow=[");
% for ii=1:length(lastrow)
%     ti=strcat(ti,num2str(lastrow(ii)),",");
% end
% ti=strcat(ti,"]")
% title(ti)






