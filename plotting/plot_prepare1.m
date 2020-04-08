tticks=[1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2];
la=1+tticks;

%sort mmarray into bins
mm2=mmarray;
mm2(:,:)=max(max(mmarray));
for i=length(la):-1:1
    la(i)
    nn=find(mmarray<la(i));
    size(nn)
    mm2(nn)=la(i);
end

%now change everything , shift by 1
mm2=mm2-1;
la=la-1;

cp=colormap(parula);

colormap([1 1 1; ...
    0.9769    0.9839    0.0805;...
    0.9628    0.9373    0.1265;...
    0.9786    0.8386    0.1766;...
    0.3406    0.8008    0.4789;...
    0.2440    0.4358    0.9988]);