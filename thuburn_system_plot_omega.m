
function thuburn_system_plot_omega()

global_space;

%plotting
LL=2.0;
MS=10;

%get Thuburn's matrices
[nonstiffA,stiffB] = thuburn_system();

nonstiffA = -nonstiffA/imm;
stiffB = -stiffB/imm;

a = nonstiffA + stiffB;

%vectors and diag matrix with eigvals
[aa,bb] = eig(a);
kk = diag(bb);
%figure(dispom);
%plot(real(kk), imag(kk), '*'); title('omega');

branch1=zeros(mlev-1,1);
branch2=zeros(mlev-1,1);
branch3=zeros(mlev-1,1);

%semi-analytic values of a^2 in eval problem p_zz+ap=0
kvals_load;
%analytic branches using kvals
for i=1:length(kvals)
    fun = @(x) dispR(x,kvals(i));
    root = -2;
    x1 = fzero(fun,root);
    branch1(i) = x1;
    
    root = -1e-3;
    x2 = fzero(fun,root);
    branch2(i) = x2;
    
    root = -1e-6;
    x3 = fzero(fun,root);
    branch3(i) = x3;
end

figure(dispfig); semilogy(kvals*D/pi,abs(branch1), 'b-','LineWidth', LL); hold on;
figure(dispfig); semilogy(kvals*D/pi,abs(branch2), 'b-','LineWidth', LL); hold on;
figure(dispfig); semilogy(kvals*D/pi,abs(branch3), 'b-','LineWidth', LL); hold on;
%title('system [w phi, u v sigma], p=0 on top')

plotom = [];
wavemode=zeros(numvar,1);
for ii=1:numvar
    
    ev = real(kk(ii));
    
    %find n from w :
    vec = imag(aa(wrange,ii));
    
    vecs = real(aa(srange,ii));
    vecpt = DiffPhi*real(aa(phirange,ii));
    
    %pressure in theta
    vecpress = pressureR_mid'.*(diag(1./sigmaR_mid)/kap1*vecs - diag(1./dphiR_mid)/kap1*vecpt);
    
    %pressure in z is given by
    %pressure|z = pressure|theta - z*(dp^ref/dz)
    %vecp_in_z = vecpress - zz.* pRz;
    
    vect_in_z = vecpress - thetam'.* ( (thetam').^(-1/kappa-1)*T^(-1/kappa)*(-1/kappa) );
    
    %vecpress = vecpress.*exp(-zz*(A+B)/2);
%%%%%%%%%% prev revision used vecpress
    %vec = vecpress;
    %vec = vect_in_z.*exp(-zz*(Ac+Bc)/2);
    nval = abs(obtainN(vec)) - 1;
    
    %assign smallest wavenumbers by eye for mlev20
    if 1
    %low wave numbers    
    if(ii==40)
        nval=1/3;
    elseif(ii==41)
        nval=1/3;
    elseif(ii==42)
        nval=1/3;
    elseif(ii==39)
        nval=1/3;
    elseif(ii==81)
        nval=1/3;
    %high wavenumbers rossby
    elseif(ii==99)
        nval=19;
    elseif(ii==93)
        nval=17;
    %high wavenumbers gw, red
    elseif(ii==77)
        nval=19;
    elseif(ii==78)
        nval=16;
    elseif(ii==79)
        nval=18;  
    elseif(ii==80)
        nval=17;
    %high wavenumbers gw, blue
    elseif(ii==73)
        nval=19;
    elseif(ii==74)
        nval=16;
    elseif(ii==75)
        nval=18;  
    elseif(ii==76)
        nval=17;    
    end
    end
    wavemode(ii)=nval;
    %plot individual eig functions

    omscaled = real(ev);
    
     if ((abs(omscaled) < 1e-3)&&(abs(omscaled) > 1e-5))
% [ii,nval]
% pause
         %      figure(10000+ii)
%      plot(vec, 'b')
%      title(strcat('nval here is ',int2str(nval)))
     end

    if(omscaled > 0)
        figure(dispfig); semilogy(abs(nval),omscaled,'r*','MarkerSize',MS,'LineWidth',LL); hold on;
    elseif(omscaled<0)
        figure(dispfig); semilogy(abs(nval),abs(omscaled),'bd','MarkerSize',MS,'LineWidth',LL); hold on;
    else
        figure(dispfig); semilogy(abs(nval),omscaled,'r.','MarkerSize',MS,'LineWidth',LL); hold on;
    end
    
end

ti=strcat('space-discr-mlev',int2str(mlev),'.mat');
save(ti,'aa','kk','nonstiffA','stiffB','wavemode');

%print('myfig','-dpdf')

set(gca,'FontSize',18)
xlabel('vertical mode $n$','Interpreter', 'latex')
ylabel('frequency $\omega$','Interpreter', 'latex')

end






