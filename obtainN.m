%counting zeros to get a wavenumber

function res = obtainN(vec)
ll = length(vec);
if abs(norm(vec))<1e-14
    res = 0;
else
    if vec(1) > 0 
        si = 1;
    else
        si =-1;
    end
    prev = si; ch = 1;
    for ii = 2:ll
    if ( vec(ii)*prev > 0 )
    else
        ch = ch+1;
        prev = sign(vec(ii));
    end
    end
    
    res = ch*si;
end
end