%global variables/functions/profiles

global T R cp g pa kappa cr2 k K2 f bbeta Nr2 Ac Bc Gamma D kap1 mlev EC resolution wavenumber;
global l rt imm;
global wrange;
global dispfig dispom;

EC=40*10^6;
D=10^4;

bbeta = 1.619e-11;
f = 1.031e-4;
g = 9.80616;
R = 287.05;
cp = 1005.0;
kappa = R/cp;
T = 250;

numvar=5*mlev;

%original Thuburns number k = 2*pi/1000000, is wavenumber 40 on unit-radius sphere
%since Earth's radius is 40*10^6
%k = 2*pi/1000000 = 2*pi*40 / 40*10^6;
%smallest resolved wavenumber for 1 degree sem,
%wavenumber 180 on unit-radius earth

%wavenumber k=2pi/T, T is wavelength
%if input wavenumber is 180, then wavelength is 2dx,
%if input = 90, wavelength is 4dx
%so T = 360/wavenumber * earthcirc/12/res,
% dx = earthcirc/12/res,
% 12 here is 3*4, 3 = np-1, np=4 from usual homme setup
% and 4 is 4 cube faces.
wavelength=360 / wavenumber * EC / 12 / resolution;
k=2*pi/wavelength;

%rest of parameters
l = 0;
K2 = k*k+l*l;
kap1 = 1-kappa;
cr2 = R*T/kap1;
Nr2 = g*g *kappa/R/T;
Ac = Nr2/g;
Bc = g/cr2;
Gamma = (Bc-Ac)/2;
rt = R*T;
pa = g/R/T;

%reference profiles and constants
thetatop = thetaR_inz(D);
thetabottom = thetaR_inz(0);
dtheta = (thetatop-thetabottom)/mlev;

thetaint = [thetabottom:dtheta:thetatop];
thetam = thetaint(1:end-1)+dtheta/2;
imm = complex(0,1);

numvar  = mlev*5;
%var = zeros(numvar,1);
nonstiffA = zeros(numvar);
stiffB = nonstiffA;

%Equation 1: -i om u = fv+ik eta/K^2 u -ik p/ pho_r - ik phi: %%%%%%%%%%%%%%%%%%%%%%%
rhoRm = densityR_intheta(thetam);
rhoRi = densityR_intheta(thetaint(2:end));

sigmaR_mid = densityR_intheta(thetam).*dpdtR_intheta(thetam)/g;
sigmaR_int = densityR_intheta(thetaint(2:end)).*dpdtR_intheta(thetaint(2:end))/g;
pressureR_mid = pressureR_intheta(thetam);
dphiR_mid = dpdtR_intheta(thetam);
dphiR_int = dpdtR_intheta(thetaint(2:end));

AverSigma = buildAverSigma(mlev);
DiffPressure = buildDiffPressure(mlev,dtheta);
AverPhi = buildAverPhi(mlev);
DiffPhi = buildDiffPhi(mlev,dtheta);

urange=1:mlev;
vrange=urange+mlev;
wrange=urange+2*mlev;
phirange=urange+3*mlev;
srange=urange+4*mlev;

%%%% for plotting omega
zz = zoft(thetam');
pRz = pressureRz_inz(zz);


imm = complex(0,1);

%good
function res = pressureR_inz(z)
global R T pa;
res = exp(-pa*z)*R*T;
end
%good
function res = pressureRz_inz(z)
global R T pa;
res = -pa*exp(-pa*z)*R*T;
end
%good
function res = densityRz_inz(z)
global pa;
res = -pa*exp(-pa*z);
end
%good
function res = densityR_inz(z)
global pa;
res = exp(-pa*z);
end
%good
function res = thetaR_inz(z)
global T kappa;
res = T*(pressureR_inz(z)).^(-kappa);
end

function res = pressureR_intheta(t)
res = pressureR_inz(zoft(t));
end
function res = pressureRz_intheta(t)
res = pressureRz_inz(zoft(t));
end
function res = densityRz_intheta(t)
res = densityRz_inz(zoft(t));
end
function res = densityR_intheta(t)
res = densityR_inz(zoft(t));
end

function res = zoft(t)
global R T g kappa;
res = R*T/g/kappa*log((R*T)^kappa*t/T);
end
function res = dpdtR_intheta(t)
global R T kappa;
res = R*T/kappa./t;
end

%computing dphi on midpoints, with zero BC on bottom
% take mlev interface values (from 2nd to mlev+1 interface)
% and obtain mvel mid level values
%verified
function res = buildDiffPhi(mlev,de)
res = zeros(mlev,mlev);
for i=2:mlev
    res(i,i-1) = -1; res(i,i) = 1;
end
res(1,1) = 1;
res = res/de;
end

%computing phi on midpoints, with zero BC on bottom
% take urange interface values and obtain mvel mid level values
% so, we dont solve w, \phi on top
%verified
function res = buildAverPhi(mlev)
res = zeros(mlev,mlev);
for i=2:mlev
    res(i,i-1) = 1/2; res(i,i) = 1/2;
end
res(1,1) = 1/2;
end

%computing dQ/deta on interfaces
% take mlev middle values of Q and obtain mlev interface values
% of dQ for 2:END interfaces
% with zero BC on top, Q=0
%verified
function res = buildDiffPressure(mlev,de)
res = zeros(mlev,mlev);
for i=1:mlev - 1
    res(i,i) = -1; res(i,i+1) = 1;
end
res(mlev,mlev) = -2;
res = res/de;
end

%computing sigma on 2:END interfaces without BC
% the value for the top interface is the same as the value for
% the closest mindpoint -- extrapolation.
%verified
function res = buildAverSigma(mlev)
res = zeros(mlev,mlev);
for i=1:mlev-1
    res(i,i) = 1/2; res(i,i+1) = 1/2;
end
res(mlev,mlev) = 1;
end






