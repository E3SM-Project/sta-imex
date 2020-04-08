function res = IMEXRKstabmat(N,S,d,locdt,A,Ahat,b,bhat,r)
M = zeros(r,d,d);
Id = eye(d,d);
M(1,1:d,1:d) = Id;
temp = zeros(d,d);
temp(:,:) = M(1,:,:);
R = Id + locdt *(b(1)*N + bhat(1)*S)*temp(:,:); %%%

for j = 2:r
    summ = eye(d,d);
    for l=1:j-1
        temp(:,:) = M(l,:,:);
        summ = summ + locdt *(A(j,l)*N + Ahat(j,l)*S)*temp(:,:); %%%
    end
    M(j,1:d,1:d) = inv(Id - locdt *Ahat(j,j)*S)*summ;  %%%
    temp(:,:) = M(j,:,:);
    R = R+ locdt * (b(j)*N + bhat(j)*S)*temp(:,:); %%%
end
res = R;
end