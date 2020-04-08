function paper10101
global A Ahat b bhat c chat r
    r = 6; coefa1 = 1/4; coefa2 = 1/6; coefa3 = 3/8; coefa4 = 1/2; coefa5 = 1;
    A = [0 0 0 0 0 0 ; ...
         coefa1 0 0 0 0 0 ; ...
         0 coefa2 0 0 0 0 ; ...
         0 0 coefa3 0 0 0 ; ...
         0 0 0 coefa4 0 0 ; ...
         0 0 0 0 coefa5 0];
    c = A*ones(r,1); b = A(r,:)';
    
    ahat1 = 0 ;
    dhat1 = coefa1 ;
    
    ahat2 = 0 ; 
    dhat2 = coefa2 ;
    
    ahat3 = 0 ; 
    dhat3 = coefa3 ;
    
    ahat4 = 0 ; 
    dhat4 = coefa4 ;

     ee=0.02;
    ahat51 = 1/2-ee;
    ahat52 = 0;
    ahat53 = 0;
    ahat54 = 0;
    ahat55 = 0;
    dhat5 =  1/2+ee;



    Ahat = [ 0 0 0 0 0 0 ; ahat1 dhat1 0 0 0 0 ; ahat2 0 dhat2 0 0 0 ; ...
             ahat3 0 0 dhat3 0 0 ; ahat4 0 0 0 dhat4 0 ; ...
             ahat51 ahat52 ahat53 ahat54 ahat55 dhat5];
    chat = Ahat*ones(r,1); bhat = Ahat(r,:)';
end