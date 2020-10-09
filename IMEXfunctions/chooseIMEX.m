function chooseIMEX(which)

%0 is some ARS
%1 is ttype7
%2 is S1S2
%3 is S3

if(which==0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEFORE PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PAPER NOTATIONS 
elseif(which==10090)
%original ttype 9      s53
    paper10090;
elseif(which==10091)
%modif ttype 9 Feb01
    paper10091;
elseif(which==10092)
%ttype 9 with last BE
    paper10092;
elseif(which==10096)
%ttype 9 experimental
    paper10096;    
elseif(which==10100)
%original ttype 10 S4 , M2a
    paper10100;
elseif(which==10101)
%modif ttype 10 Feb 01, M2b
    paper10101;
elseif(which==10102)
%modif ttype 10 AS version, M2c
    paper10102;
elseif(which==10103)
%ttype 10 with last BE, M2be
    paper10103;
elseif(which==10104)
%ttype 10 with last CN, M2cn
    paper10104;
elseif(which==10105)
%ttype 10 with last CN+offcenter, M2cno
    paper10105; 
elseif(which==10106)
%experimental
    paper10106;    
    
elseif(which==10232)
%ARK2(2,3,2) as in Rokhzadi2018
    paper10232;
elseif(which==102321)
%ARK2(2,3,2) different version, with a32=0.85
    paper102321;
elseif(which==10532)
%developed scheme in Rokhzadi2018 imex-ssp2(2,3,2)
    paper10532;

elseif(which==10030)
%S3 scheme on webpage, M1 scheme in paper
    paper10030;
elseif(which==10060)
%S6 scheme on webpage
    paper10060;

end

end