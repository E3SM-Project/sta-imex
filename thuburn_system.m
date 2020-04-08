% Produces N and S matrices in Thuburn's system. For omega, do A=N+S,
% for imex stability, nonstiffA=N is explicit/nonstiff,
% stiffB=S is implicit/stiff part.

% Resolution is given in NE notations, that is as in NE in homme.
% Hor. wavenumber k is computed from resolution, wavenumber input, and
% assumption that np=4, that is, 4 GLLs per element. Wavenumber input can
% be anything, but most interesting/challenging is 180 (corresponding to
% 2dx wavelength) and 90 (corresponfding to 4dx wavelength).

% Domain usually should be 10K as equivalent height, as in Thuburn.
% EC is Earth's circumference for k, EC = 40*10^6

%%%%%%%%%%% Lagrangian system
%%% Staggering [w \phi, u v sigma]
%%% with BC p=0 on top.
%%% Here sigma term in u eqn for p/p^r term is not interpolated to interfaces
%%% and back to midlevels. Instead just midlevel sigma values are used, thus,
%%% omitting BC condition for sigma on the top.
%%% Also, in w eqn sigma term is also treated without BC,
%%% but dsigma term in dp is treated with matrix Ap (which implies p=0 on top).
%%% The only difference betweed this script and the one with _dphi_bc
%%% is in this line in w eqn: <-- matrix Ap instead of matrix NS
%%%             - diag(1./sigmaR_int)*matrNS*diag(pressureR_mid./sigmaR_mid)/kap1;

%%% not Thuburn's result

function [nonstiffA,stiffB] = thuburn_system()

global_space;

ilev=mlev+1;
slev=mlev-1;

%ordering of variables in A, B matrices
% u v w \phi sigma

for ii = urange
    %term fv  VERIFIED
    nonstiffA(ii, ii+mlev) = nonstiffA(ii, ii+mlev) + f;
    %term ik beta /K^2 u VERIFIED
    nonstiffA(ii,ii) = nonstiffA(ii,ii) + imm*k*bbeta/K2;
end
%-ik \phi term
nonstiffA(urange,phirange) = ...
    nonstiffA(urange,phirange) - imm*k*AverPhi;
%part of term -ik p/rho^r, based on phi
nonstiffA(urange,phirange) = ...
    nonstiffA(urange,phirange)    + imm*k/kap1*diag(rt./dphiR_mid)*DiffPhi;
%sigma
nonstiffA(urange,srange) = ...
    nonstiffA(urange,srange) - imm*k*rt/kap1*diag(1./sigmaR_mid);

%Equation 2: - i om v = -fu+ i k beta/k^2 v - il*p/pho_r - il phi %%%%%%%%%%%%
for ii = vrange
    %term -fu VERIFIED
    nonstiffA(ii, ii-mlev) = nonstiffA(ii, ii-mlev) - f;
    %term ik beta /K^2 v VERIFIED
    nonstiffA(ii,ii) = nonstiffA(ii,ii) + imm*k*bbeta/K2;
end

%third equation:-i om w =...
%Let's start collecting terms for eqn 3: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%term -g*sigma/sigmaR
stiffB(wrange, srange) = ...
    stiffB(wrange, srange) - g*diag(1./sigmaR_int)*AverSigma;

%pressure term based on sigma
stiffB(wrange,srange) = ...
    stiffB(wrange,srange)  ...
    - diag(1./sigmaR_int)*DiffPressure*diag(pressureR_mid./sigmaR_mid)/kap1;

%pressure term based on phi
stiffB(wrange,phirange) = ...
    stiffB(wrange,phirange)  ...
    + diag(1./sigmaR_int)*DiffPressure*diag(pressureR_mid./dphiR_mid)*DiffPhi/kap1 ;

%equation 4: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -i om phi =  g w - etadot rhi^r_\eta
for ii = urange
    phiind = ii + 3*mlev;
    %shift to w index
    wind = 2*mlev+ii;
    stiffB(phiind, wind) = ...
        stiffB(phiind, wind) + g;
end

%equation 5: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=urange
    aind = 4*mlev+ii;
    uind = ii;
    nonstiffA(aind,uind) = ...
        nonstiffA(aind,uind) - imm*k*sigmaR_mid(ii);
end


end %end of thuburn_system


