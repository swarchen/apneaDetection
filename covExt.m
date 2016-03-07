function [covTA,cov]=covExt(par,LEN,STAGE,statePSG,THO,ABD,sp,ep,CFlow)
PAR=par.par;
CW=par.cov;
SLIDE=par.slide;
NOR=par.st.nor;
OSA=par.st.osa;
CSA=par.st.csa;
MSA=par.st.msa;
HYP=par.st.hyp;
fs=par.fsDown;
qt=par.qt;
stdR=1;%normal ratio
curR=1;%current ratio
preR=1;
unk=zeros(length(statePSG),1);
nor=zeros(length(statePSG),1);
evn=zeros(length(statePSG),1);
covTA=zeros(length(statePSG),1);

for ss=sp:ep%PAR+1:LEN-CW+1%fix(CW/2)+1:LEN-fix(CW/2)%fix(FR/2)+1:LEN-fix(FR/2)
    stg=STAGE(fix(ss/SLIDE));
    cwR=fs*(ss-1)/SLIDE+1:fs*(ss+CW-1)/SLIDE;
    cwTHO=THO(cwR);
    cwABD=ABD(cwR);
    cwCFL=CFlow(cwR);
    if stg==11
        covTA(ss,1)=nan;
        unk(ss,1)=0;
        nor(ss,1)=0;
        evn(ss,1)=0;
    else
        isosa=( find( statePSG(ss:ss+CW-1)==OSA ) );
        iscsa=( find( statePSG(ss:ss+CW-1)==CSA ) );
        ismsa=( find( statePSG(ss:ss+CW-1)==MSA ) );
        ishyp=( find( statePSG(ss:ss+CW-1)==HYP ) );
        isnor=( find( statePSG(ss:ss+CW-1)==NOR ) );
        covTA(ss,1)=dot(cwTHO,cwABD)/norm(cwTHO)/norm(cwABD);
        if( length(isosa)==CW )
            unk(ss,1)=0;
            nor(ss,1)=0;
            evn(ss,1)=OSA;
        elseif( length(iscsa)==CW )
            unk(ss,1)=0;
            nor(ss,1)=0;
            evn(ss,1)=CSA;
        elseif( length(ishyp)==CW )
            unk(ss,1)=0;
            nor(ss,1)=0;
            evn(ss,1)=HYP;
        elseif( length(ismsa)==CW )
            unk(ss,1)=0;
            nor(ss,1)=0;
            evn(ss,1)=MSA;
        elseif( length(isnor)==CW )
            unk(ss,1)=0;
            nor(ss,1)=1;
            evn(ss,1)=0;
        else
            unk(ss,1)=1;
            nor(ss,1)=0;
            evn(ss,1)=0;
        end
    end
end
cov.unk=unk;
cov.nor=nor;
cov.evn=evn;
covTA(1:sp-1)=nan;
covTA(ep+1:end)=nan;
end
