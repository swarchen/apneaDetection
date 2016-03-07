function [mfTHO,mfABD,rfTHO,rfABD,fq]=FrExt(par,LEN,STAGE,statePSG,THO,ABD,sp,ep)
FR=par.fr;
SLIDE=par.slide;
NOR=par.st.nor;
OSA=par.st.osa;
CSA=par.st.csa;
MSA=par.st.msa;
HYP=par.st.hyp;
qt=par.qt;
stdR=1;%normal ratio
curR=1;%current ratio
preR=1;
unk=zeros(length(statePSG),1);
nor=zeros(length(statePSG),1);
evn=zeros(length(statePSG),1);
mfTHO=zeros(length(statePSG),1);
mfABD=zeros(length(statePSG),1);
rfTHO=zeros(length(statePSG),1);
rfABD=zeros(length(statePSG),1);

for ss=sp:ep
    stg=STAGE(fix(ss/SLIDE));
    fqR=par.fsDown*(ss-1)/SLIDE+1:par.fsDown*(ss+FR-1)/SLIDE;%fs/SLIDE*(ss-fix(FR/2))+1:fs/SLIDE*(ss+fix(FR/2));
    fqTHO=THO(fqR);
    fqABD=ABD(fqR);
    fqTHO=fqTHO-mean(fqTHO);
    fqABD=fqABD-mean(fqABD);
    downfac=20;
    fs=par.fsDown;%par.fs/downfac;
    fqTHO=fqTHO'-smooth(fqTHO,5*fs,'loess');
    fqABD=fqABD'-smooth(fqABD,5*fs,'loess');
    f1 = 0; 
    f2 = 4;     % in hertz
    m = 128;
    df=fs/m;
    f_axis=df:df:fs;
    fTHO=abs(fft(fqTHO,m)).^2;%(abs(czt(fqTHO,m,w,s))).^2;
    fABD=abs(fft(fqABD,m)).^2;%(abs(czt(fqABD,m,w,s))).^2;
    Erg=fix(0.8/fs*m):fix(1.5/fs*m);
    Erg1=fix(0.1/fs*m):fix(0.8/fs*m); 
    if stg==11
        mfTHO(ss,1)=nan;
        mfABD(ss,1)=nan;
        rfTHO(ss,1)=nan;
        rfABD(ss,1)=nan;
        unk(ss,1)=0;
        nor(ss,1)=0;
        evn(ss,1)=0;
    else
        isosa=( find( statePSG(ss:ss+FR-1)==OSA ) );
        iscsa=( find( statePSG(ss:ss+FR-1)==CSA ) );
        ismsa=( find( statePSG(ss:ss+FR-1)==MSA ) );
        ishyp=( find( statePSG(ss:ss+FR-1)==HYP ) );
        isnor=( find( statePSG(ss:ss+FR-1)==NOR ) );
        mft=(find(fTHO(1:fix(m/2))==max(fTHO(1:fix(m/2)))));
        mfa=(find(fABD(1:fix(m/2))==max(fABD(1:fix(m/2)))));
        mfTHO(ss,1)=(mft)*df;
        mfABD(ss,1)=(mfa)*df;
        rfTHO(ss,1)=sum(fTHO(Erg))/sum(fTHO(Erg1));
        rfABD(ss,1)=sum(fABD(Erg))/sum(fABD(Erg1));
        if( length(isosa)==FR )            
            unk(ss,1)=0;
            nor(ss,1)=0;
            evn(ss,1)=OSA;
        elseif( length(iscsa)==FR )            
            unk(ss,1)=0;
            nor(ss,1)=0;
            evn(ss,1)=CSA;
        elseif( length(ishyp)==FR )
            unk(ss,1)=0;
            nor(ss,1)=0;
            evn(ss,1)=HYP;
        elseif( length(ismsa)==FR )
            unk(ss,1)=0;
            nor(ss,1)=0;
            evn(ss,1)=MSA;
        elseif( length(isnor)==FR )
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
fq.unk=unk;
fq.nor=nor;
fq.evn=evn;
mfTHO(1:sp-1)=nan;
mfABD(1:sp-1)=nan;
mfTHO(ep+1:end)=nan;
mfABD(ep+1:end)=nan;
rfTHO(1:sp-1)=nan;
rfABD(1:sp-1)=nan;
rfTHO(ep+1:end)=nan;
rfABD(ep+1:end)=nan;
end
