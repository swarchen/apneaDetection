function [rtoTHO,rtoABD,rt]=rtExt(par,LEN,STAGE,statePSG,THO,ABD,sp,ep)
PAR=par.par;
CAR=par.car;
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
rtoTHO=zeros(length(statePSG),1);
rtoABD=zeros(length(statePSG),1);

for ss=sp:ep
    stg=STAGE(fix(ss/SLIDE));
    preR=fs*(ss-PAR-1)/SLIDE+1:fs*(ss-1)/SLIDE;
    if( ((statePSG(ss)==OSA)||(statePSG(ss)==CSA)||(statePSG(ss)==HYP))&&(statePSG(ss-1)==NOR) )
        stdR=fs*(ss-PAR-1)/SLIDE+1:fs*(ss-1)/SLIDE;
    end
    curR=fs*(ss-1)/SLIDE+1:fs*(ss+CAR-1)/SLIDE;
    curTHO=THO(curR);
    curABD=ABD(curR);
    stdTHO=THO(stdR);
    stdABD=ABD(stdR);
    preTHO=THO(preR);
    preABD=ABD(preR);
    if stg==11
        rtoTHO(ss,1)=nan;
        rtoABD(ss,1)=nan;
        unk(ss,1)=0;
        nor(ss,1)=0;
        evn(ss,1)=0;
    else
        isosa=( find( statePSG(ss:ss+CAR-1)==OSA ) );
        iscsa=( find( statePSG(ss:ss+CAR-1)==CSA ) );
        ismsa=( find( statePSG(ss:ss+CAR-1)==MSA ) );
        ishyp=( find( statePSG(ss:ss+CAR-1)==HYP ) );
        isnor=( find( statePSG(ss:ss+CAR-1)==NOR ) );
        if( length(isosa)==CAR )
            rtoTHO(ss,1)=quantile(abs((curTHO)),qt)/quantile(abs((stdTHO)),qt);%mean(abs(curTHO))/mean(abs(stdTHO));
            rtoABD(ss,1)=quantile(abs((curABD)),qt)/quantile(abs((stdABD)),qt);%mean(abs(curABD))/mean(abs(stdABD));
            unk(ss,1)=0;
            nor(ss,1)=0;
            evn(ss,1)=OSA;
        elseif( length(iscsa)==CAR )
            rtoTHO(ss,1)=quantile(abs((curTHO)),qt)/quantile(abs((stdTHO)),qt);%mean(abs(curTHO))/mean(abs(stdTHO));
            rtoABD(ss,1)=quantile(abs((curABD)),qt)/quantile(abs((stdABD)),qt);%mean(abs(curABD))/mean(abs(stdABD));
            unk(ss,1)=0;
            nor(ss,1)=0;
            evn(ss,1)=CSA;
        elseif( length(ishyp)==CAR )
            rtoTHO(ss,1)=quantile(abs((curTHO)),qt)/quantile(abs((stdTHO)),qt);%mean(abs(curTHO))/mean(abs(stdTHO));
            rtoABD(ss,1)=quantile(abs((curABD)),qt)/quantile(abs((stdABD)),qt);%mean(abs(curABD))/mean(abs(stdABD));
            unk(ss,1)=0;
            nor(ss,1)=0;
            evn(ss,1)=HYP;
        elseif( length(ismsa)==CAR )
            rtoTHO(ss,1)=quantile(abs((curTHO)),qt)/quantile(abs((stdTHO)),qt);%mean(abs(curTHO))/mean(abs(preTHO));
            rtoABD(ss,1)=quantile(abs((curABD)),qt)/quantile(abs((stdABD)),qt);%mean(abs(curABD))/mean(abs(preABD));
            unk(ss,1)=0;
            nor(ss,1)=0;
            evn(ss,1)=MSA;
        elseif(length(isnor)==CAR)
            rtoTHO(ss,1)=quantile(abs((curTHO)),qt)/quantile(abs((preTHO)),qt);%mean(abs(curTHO))/mean(abs(preTHO));
            rtoABD(ss,1)=quantile(abs((curABD)),qt)/quantile(abs((preABD)),qt);%mean(abs(curABD))/mean(abs(preABD));
            unk(ss,1)=0;
            nor(ss,1)=1;
            evn(ss,1)=0;
        else
            rtoTHO(ss,1)=quantile(abs((curTHO)),qt)/quantile(abs((preTHO)),qt);%mean(abs(curTHO))/mean(abs(preTHO));
            rtoABD(ss,1)=quantile(abs((curABD)),qt)/quantile(abs((preABD)),qt);%mean(abs(curABD))/mean(abs(preABD));
            unk(ss,1)=1;
            nor(ss,1)=0;
            evn(ss,1)=0;
        end
    end
end
rt.unk=unk;
rt.nor=nor;
rt.evn=evn;
rtoTHO(1:sp-1)=nan;
rtoABD(1:sp-1)=nan;
rtoTHO(ep+1:end)=nan;
rtoABD(ep+1:end)=nan;
end