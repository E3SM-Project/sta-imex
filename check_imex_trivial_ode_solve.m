%check imex timestepping for oaa 
%analyt function is
%f(t) = t^3 + sin^2(20*t)

%a is nonstiff, b is stiff
clear all
numvar=1;
global nn ss
nn=0.1; ss=0.9;
t0 = 0; 
t1 = 1;
dt=[1 0.1 0.01 0.001 ];%0.0001 0.00001 0.000001];


%%% CHOOSE IMEX %%%%%%%%%%%%%%%%%%%%%%%%%%

global A Ahat b bhat c chat gamma r

%0 is some ARS
%1 is ttype7
%2 is S1S2
%3 is S3
which=10102;
chooseIMEX(which);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%1 step
for i=1:length(dt)
    dtl=dt(i);
    step = IMEXRKstabmat(nn,ss,1,dtl,A,Ahat,b,bhat,r);
    %[fa(t0+dtl), (fa(t0)* step), abs(fa(t0+dtl) - (fa(t0)* step) )]
    %[abs(fa(t0+dtl) - (fa(t0)* step) )]
    %%% for CN only
    %cn = (1+dtl/2*ss)/(1-dtl/2*ss);
    %[cn - step]
end
aa=1;

%over the same interval [1 11]
for i=1:length(dt)
    dtl=dt(i);
    
    initf = fa(t0);
    currf = initf;
    nsteps = (t1-t0)/dtl;
    for j=1:nsteps
       step = IMEXRKstabmat(nn,ss,1,dtl,A,Ahat,b,bhat,r);
       currf = currf * step;
    end
    %[fa(t1), currf , fa(t1) - currf]
    [fa(t1) - currf]
end


function res=fa(t)
global nn ss
res = exp((nn+ss)*t);
end





