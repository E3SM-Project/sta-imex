function paper10232
global A Ahat b bhat c chat r

r = 3;

ssq = sqrt(2);
gamma = 1 - 1/ssq;
delta = 3 + 2*ssq;

A = [0 0 0 ; ...
    2*gamma 0 0 ; ...
    1-delta/6 delta/6 0];

b = [1/2/ssq 1/2/ssq gamma]';
c = A*ones(r,1);

Ahat = [0 0 0 ;...
    gamma gamma 0 ;...
    1/2/ssq 1/2/ssq gamma];

bhat = Ahat(r,:)';
chat = Ahat*ones(r,1);
    
end