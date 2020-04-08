%based on analytical eigenvalues \lam of the problem
%p_zz + \lam p = 0. p=0 on top (z=D), p_z+Bp=0 on bottom (z=0), 
%where \lam is found
%from an exact eqn for \lam=k^2, \tan(kD)=k/B numerically,
%here we find 3 branches of acoustic, gravity, and rossby waves,
%similarly to Thuburn's result with \omegas.

function main

clear all;

global resolution wavenumber mlev dispfig dispom;

%res and wavenumber are dummies here
resolution=1;
wavenumber=360*40/12;

%mlev is meaningful
mlev=20;

global_space;


% first, find first 20 values for sqrt(a)
% for p_zz+ap=0, p_z+(B/2-A/2)p=0 at z=0,D
ta = zeros(mlev,1);
fun2 = @dummy;
% solutions will obey tan(sqrt(a)*D)=sqrt(a)/(B/2-A/2)
% scale out D: solve instead
% tan(sqrt(ta))=sqrt(ta)/(B/2-A/2)/D
% 1/(B/2-A/2)/D ~= 3.4, so, very steep slope
% also, ta=0 is not a solutions satisfying BC
% for large values, ta is very close to (2n+1)*pi/2
for ii=1:mlev
    if (ii == 1)
        xzero = [0.05,pi/2*.95];
    else
        xzero = [pi/2*.5, pi/2*0.9999999] + (ii-1)*pi;
    end
    ta(ii) = fzero(fun2,xzero);
end

%these are sqrt(a) from p_{zz} + a*p=0 .
kvals = ta/D;
%let's not compute coefficient a... leave sqrt(a) here instead.
%kvals^2 = C(omega)-(B-A)^2/4, so, let's find 5 omegas:

kvals

%ro = rossby
%gn = gravity
%an = acoustics
ro = zeros(mlev,1);
gn = ro; an = ro;
for i=1:20
    ava = kvals(i);
    
    fun3 = @(x) dispR(x,ava);
    
    gn(i)=fzero(fun3,-1e-3);
    an(i)=fzero(fun3,[-2.3, -0.01]);
    ro(i)=fzero(fun3,-1e-6);
    
end

figure(23)
title("Three wave branches")
semilogy(kvals,abs(ro),'b*-'); hold on;
semilogy(kvals,abs(an),'r*-'); hold on;
semilogy(kvals,abs(gn),'black*-'); hold on;

%print('hommesp','-dpdf')
% ro
% an
% gn

end

function res = dummy(k)
global Ac Bc D;
ab=(Bc-Ac)/2;
res = tan(k) - k/ab/D;
end

function res = C(om)
global cr2 Nr2 K2 k bbeta f;
res = ( 1/cr2 + (K2+k*bbeta./om)./(f*f-(om+k*bbeta/K2).^2) ).*(om.*om-Nr2);
end


