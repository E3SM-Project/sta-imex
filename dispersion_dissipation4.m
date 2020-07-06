%clear all;

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
%which=10105;

DT=50;

global_space;
global_space_imex;
[nonstiffA,stiffB] = thuburn_system();
        
chooseIMEX(which);

%spacetime op
ime = IMEXRKstabmat(nonstiffA,stiffB,numvar,DT,A,Ahat,b,bhat,r);
[imexvector,imexeig]=eig(ime);
imexeigval=diag(imexeig);

%ti=strcat('for-spacetime-dispersion/method-mlev',int2str(mlev),'-',int2str(which),'-dt',num2str(DT),'.mat');
%save(ti,'imexA','imexAvec','imexAeig','DT');

%plotting
LL=2.0;
MS=10;

%eigenvectors, eigenvalues of space discr, nonstiffA, stiffB
ti=strcat('space-discr-mlev',int2str(mlev),'.mat');
aafile=load(ti);

figg=22400000;

costamb=1;

%total matr A and parts N,S for A=N+S for spacetime op
imexfile=load(ti);

%eigenvectors of space discr
oma=aafile.aa;
%eigenvals of space discr
speig=aafile.kk;
%nonstiff matrix
nonstiff=aafile.nonstiffA;
%stiff matrix
stiff=aafile.stiffB;
nvals=eig(nonstiff);
svals=eig(stiff);

ss=mlev*5;
tol=1e-1;
imm=complex(0,1);


%3 types of waves, 1=acoustics, 0=gravity, 2=rossby
%define by frequencies
kkind=zeros(ss,1);
for ii=1:ss
    %rossby
    if (abs(speig(ii))<1e-5)
        kkind(ii)=2;
    elseif(abs(speig(ii))>1e-2)
        kkind(ii)=1;
    end
end

%sanity check count of waves
size(find(kkind==1))
size(find(kkind==2))

%plot eigvals of N + S and N and S
% figure(223)
% semilogy(zeros(ss,1),abs(real(speig)),'b*'); hold on;
% semilogy(zeros(ss,1),abs(real(nvals)),'kd'); hold on;
% semilogy(zeros(ss,1),abs(real(svals)),'ro');

BIG=1e10;

costmat=zeros(ss);
for indspace=1:ss
    
    %space vector
    if(~costamb)
    ee=real(oma(:,indspace));
    else
    ee=oma(:,indspace);
    end
    for ii=1:ss
        %sp time vector
        if(~costamb)
        ff=real(imexvector(:,ii));
        else
        ff=imexvector(:,ii);
        end
        %cost, done before with real parts only
        if(~costamb)
            %matlabs thing -- for arg=1.000000000000000 
            %acos(arg)=0.000000000000000e+00 + 2.107342425544702e-08i
        costmat(indspace,ii) = real(acos( dot(ee,ff)/norm(ee)/norm(ff)));
        else
        costmat(indspace,ii) = -abs( dot(ee,ff) )/norm(ee)/norm(ff);
        end
    end
end

[assign,cost]=munkres(costmat);
matchfound=assign*[1:ss]';

for indspace=1:ss
    indsptime=matchfound(indspace);
    timekind(indsptime)=kkind(indspace);
end
%sanity check count of waves
%gravity
size(find(timekind==0))
%rossby
size(find(timekind==2))

%plot matched pairs on top of each other
if 0
    for whi=1:5
        ee=real(oma(:,whi));
        v2=matchfound(whi);
        
        whichvv=imexaa(:,v2);
        figure(whi)
        %plot space vector ee and its match in spacetime, whi
        ti=strcat('space vector ',int2str(whi),...
            'kind=',int2str(kkind(whi)),...
            ' and sptime vector',int2str(v2),...
            ' and var is ',costmat(whi,v2));
        plot([1:ss],real(whichvv),'r*-',[1:ss],real(ee),'bo-'); hold on;
        plot([1:ss],real(ee.*exp(imm*speig(whi)*20)),'green*-')
        title(ti);
    end
end


for indspace=1:ss
    indsptime=matchfound(indspace);
    timekind(indsptime)=kkind(indspace);
end
%sanity check count of waves
size(find(timekind==0))
size(find(timekind==2))

if 0
%%%%%%%%%%%%% UNIT CIRCLE
figure(10000000+which)
ph=[0:0.05:2*pi];
plot(cos(ph),sin(ph),'b.','MarkerSize',MS); hold on;
%plot(real(imexbb),imag(imexbb),'bo')
set(gca,'FontSize',18)
%plot acou
ff=find(timekind==1);
plot(real(imexeigval(ff)),imag(imexeigval(ff)),'rd','MarkerSize',MS,'LineWidth',LL); hold on
%plot grav
ff=find(timekind==0);
plot(real(imexeigval(ff)),imag(imexeigval(ff)),'b*','MarkerSize',MS,'LineWidth',LL); hold on
%plot rossby
ff=find(timekind==2);
plot(real(imexeigval(ff)),imag(imexeigval(ff)),'green*','MarkerSize',MS,'LineWidth',LL); hold on
axis square
end

%%%%%%%%%%%%%% DISPERSION/DISSIPATION
imexfreq=zeros(ss,1);
imexampfactor=imexfreq;
X=ones(ss,1);

%true exp(imm*deltat*frequency) vs numerical
truefactor=exp(imm*DT*speig);

if 0
figure(55);
%acoustics
whichkind=1;
ff=find(kkind==whichkind);
plot(real(truefactor(ff)),imag(truefactor(ff)),'rd'); hold on;
ff=find(timekind==whichkind);
plot(real(imexeigval(ff)),imag(imexeigval(ff)),'r*'); hold on;

whichkind

%figure(56);
%gr
whichkind=0;
ff=find(kkind==whichkind);
plot(real(truefactor(ff)),imag(truefactor(ff)),'bd'); hold on;
ff=find(timekind==whichkind);
plot(real(imexeigval(ff)),imag(imexeigval(ff)),'b*'); hold on;

whichkind

%figure(57);
%rossby
whichkind=2;
ff=find(kkind==whichkind);
plot(real(truefactor(ff)),imag(truefactor(ff)),'gd'); hold on;
ff=find(timekind==whichkind);
plot(real(imexeigval(ff)),imag(imexeigval(ff)),'g*'); hold on;

whichkind
end



%plot frequencies
%original wavemodes from space operator
wm=aafile.wavemode;
figure(figg+which);
subplot(2,1,1);
%cycle thru space operator members
for indspace=1:ss
    %match in spacetime eigenvectors
    indsptime=matchfound(indspace);

    %exp(i*om*dt), eigval of spacetime operator
    lfactor=imexeigval(indsptime);
    
    spacevec=oma(:,indspace);
    timevec=imexvector(:,indsptime);
    
    
    
%     %are factors 1?
    newmag=abs(lfactor);
    imexampfactor(indsptime)=newmag;
%     %now phase
    %imexfreq(who)=atan(imag(lfactor)/real(lfactor))/DT;
    
    imexfreq(indsptime)=abs(real(atan(imag(lfactor)/real(lfactor))/DT));
    
    nval=wm(indspace);
    
    %here we ignore east/west propagation
    
    %rossby
    if(kkind(indspace)==2)
    semilogy(nval,real(abs(speig(indspace))),'blackd',nval,imexfreq(indsptime),'blackx',...
        'MarkerSize',MS,'LineWidth',LL);hold on;
    end
    %gravity
    if(kkind(indspace)==0)
    semilogy(nval,real(abs(speig(indspace))),'bs',nval,imexfreq(indsptime),'b+',...
        'MarkerSize',MS,'LineWidth',LL);hold on;
    end
    %acou
    if(kkind(indspace)==1)
    semilogy(nval,real(abs(speig(indspace))),'rd',nval,imexfreq(indsptime),'r*',...
        'MarkerSize',MS,'LineWidth',LL);hold on;
    end
end
%title('dispersion')
set(gca,'FontSize',18)
axis([0 mlev min(abs(speig))*0.7 max(abs(speig))*1.4])
xlabel('vertical mode $n$','Interpreter', 'latex')
ylabel('frequency $\omega$','Interpreter', 'latex')


%plot amplitudes
subplot(2,1,2)
for indspace=1:ss
    indsptime=matchfound(indspace);   
    nval=wm(indspace);
    if(kkind(indspace)==2)
    semilogy(nval,imexampfactor(indsptime),'blacko',...
        'MarkerSize',MS,'LineWidth',LL);hold on;
    elseif(kkind(indspace)==0)
    semilogy(nval,imexampfactor(indsptime),'b+',...
        'MarkerSize',MS,'LineWidth',LL);hold on;
    elseif(kkind(indspace)==1)
    semilogy(nval,imexampfactor(indsptime),'r*',...
        'MarkerSize',MS,'LineWidth',LL);hold on;
    end
end
%title('dissipation')
set(gca,'FontSize',18)
%axis([0 mlev min(imexampfactor)*0.9 1.1])
axis([0 mlev 0.01 1.1])
xlabel('vertical mode $n$','Interpreter', 'latex')
ylabel('amplification factor','Interpreter', 'latex')

x0=300;
y0=300;
width=550;
height=600
set(gcf,'position',[x0,y0,width,height])

ggg=1;
