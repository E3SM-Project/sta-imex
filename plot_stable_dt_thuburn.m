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

figure(400000);
%plotting
LL=2.0;
MS=10;

% ~ 1 km 3072
resolution=30; %high res 3000,
DTset= [250:4:600] * 30/resolution/1.5; %ne30 set, scale it
tt=length(DTset);

tol=1e-12;
%which IMEX scheme

global_space;
global_space_imex;


% whichset=[10090 10091 10100 10101 10102];
% whichlabel={"ttype9","ttype9mod","ttype10","ttype10mod","ttype10as"};
% whichcurve=["b-","black-","green-","r-","m-"];

whichset=[10100];
whichlabel={"ttype10"};
whichcurve=["g-"];


if((length(whichset) ~= length(whichlabel)) || (length(whichset) ~= length(whichcurve)))
    stop
end

ww=length(whichset);
mmarray=ones(tt,ww);

[nonstiffA,stiffB] = thuburn_system();

for jj=1:ww
    
    which = whichset(jj);
    
    chooseIMEX(which);
    
    for ii=1:tt
        
        DT=DTset(ii);
        %zb=zeros(72*5);
        imexA = IMEXRKstabmat(nonstiffA,stiffB,numvar,DT,A,Ahat,b,bhat,r);
        %imexA = IMEXRKstabmat(zb,stiffB,numvar,DT,A,Ahat,b,bhat,r);
        
        kk=eig(imexA);
        
        mmax=max(abs(kk));
        mmarray(ii,jj)=mmax;
        %define all next dt unstable if this one is unstable
%         if(mmax > 1+tol)
%             mmarray(ii:tt,jj) = 2;
%             break
%         end
        
    end
    
    %plot(DTset,mmarray(:,jj),whichcurve(jj),'DisplayName',whichlabel(jj)); hold on;
    semilogy(DTset,mmarray(:,jj)-1,whichcurve(jj),'LineWidth',LL); hold on;
    
    which
end
set(gca,'FontSize',18)

%axis([DTset(1) DTset(end) 0.9 10])

semilogy(DTset,ones(tt,1)*(1+tol-1),'black.-')

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






