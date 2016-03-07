function [state,AlgoEvNum,feat]=OvNightAnl_v14(par,doSST,LEN,rcTHO,rcABD,THO,ABD,CFlow,STAGE,statePSG,nev,svmStructNSA,svmStructOC,numTap)
P_NUM=par.name;
CAR=par.car;%/SLIDE=sec
PAR=par.par;%/SLIDE=sec
ix=1;iy=1;
sstTvlp=0;
sstAvlp=0;
SLIDE=par.slide;
fs=par.fsDown;
NOR=par.st.nor;
OSA=par.st.osa;
CSA=par.st.csa;
MSA=par.st.msa;
HYP=par.st.hyp;
data.sec=LEN/SLIDE;
data.len=LEN*par.fsDown/SLIDE;
%% SVM Factors
svNSA = svmStructNSA.SupportVectors;
alphaHatNSA = svmStructNSA.Alpha;
biasNSA = svmStructNSA.Bias;
kfunNSA = svmStructNSA.KernelFunction;
kfunargsNSA = svmStructNSA.KernelFunctionArgs;
svOC = svmStructOC.SupportVectors;
alphaHatOC = svmStructOC.Alpha;
biasOC = svmStructOC.Bias;
kfunOC = svmStructOC.KernelFunction;
kfunargsOC = svmStructOC.KernelFunctionArgs;
%%Trim Data
rcABD=rcABD(1:data.len);
rcTHO=rcTHO(1:data.len);
ABD=ABD(1:data.len);
THO=THO(1:data.len);
CFlow=CFlow(1:data.len);

num_win=(data.sec-par.pre)/SLIDE;
stdR=1;%normal ratio
curR=1;%current ratio
rtoTHO=zeros(numTap,LEN);
rtoABD=zeros(numTap,LEN);
rfTHO=zeros(numTap,LEN);
rfABD=zeros(numTap,LEN);
covTA=zeros(numTap,LEN);
fNSA=ones(numTap,1);
fOC=ones(numTap,1);

saF=0;
tF=0;
AlgoEvNum.osa=0;
AlgoEvNum.csa=0;
AlgoEvNum.hyp=0;
state=zeros(LEN,1);
for ss=PAR+1:LEN-CAR+1-(+numTap-1)
    %% Ratio Ext
    if saF==1%state(ss-1)==OSA
        stdR=par.fsDown*((ss-1)-PAR-1)/SLIDE+1:par.fsDown*((ss-1)-1)/SLIDE;
    end
    preR=par.fsDown*(ss-PAR-1)/SLIDE+1:par.fsDown*(ss-1)/SLIDE;
    for ii=0:numTap-1
        curR=par.fsDown*((ss+ii-1)-1)/SLIDE+1:par.fsDown*((ss+ii-1)+CAR-1)/SLIDE;
        curTHO(1:CAR*par.fsDown/SLIDE,ii+1)=rcTHO(1,curR)';
        curABD(1:CAR*par.fsDown/SLIDE,ii+1)=rcABD(1,curR)';
        fqTHO(1:CAR*par.fsDown/SLIDE,ii+1)=THO(1,curR)';
        fqABD(1:CAR*par.fsDown/SLIDE,ii+1)=ABD(1,curR)';
    end
    stdTHO=rcTHO(1,stdR)';
    stdABD=rcABD(1,stdR)';
    preTHO=rcTHO(1,preR)';
    preABD=rcABD(1,preR)';
    
    %% Domninant Freq & Ratio Ext
    fqTHO=fqTHO-repmat(mean(fqTHO),CAR*par.fsDown/SLIDE,1);
    fqABD=fqABD-repmat(mean(fqABD),CAR*par.fsDown/SLIDE,1);
    fs=par.fsDown;%par.fs/downfac;
    m = 64;
    Erg=fix(0.8/fs*m):fix(1.5/fs*m);
    Erg1=fix(0.1/fs*m):fix(0.8/fs*m);
    df=fs/m;
    f_axis=df:df:fs;
    rfTHO(1:numTap-1,ss)=rfTHO(2:numTap,ss-1);
    rfABD(1:numTap-1,ss)=rfABD(2:numTap,ss-1);
    fTHO=fqTHO(:,numTap);
    fABD=fqABD(:,numTap);
    fTHO=fTHO-smooth(fTHO,5*fs,'loess');
    fABD=fABD-smooth(fABD,5*fs,'loess');
    fTHO=abs(fft(fTHO,m)).^2;%(abs(czt(fqTHO,m,w,s))).^2;
    fABD=abs(fft(fABD,m)).^2;%(abs(czt(fqABD,m,w,s))).^2;
    rfTHO(numTap,ss)=log10( sum(fTHO(Erg,1))/sum(fTHO(Erg1,1)) );
    rfABD(numTap,ss)=log10( sum(fABD(Erg,1))/sum(fABD(Erg1,1)) );
    %% Covariance Ext
    cwTHO=curTHO;%THO(curR);
    cwABD=curABD;%ABD(curR);
    cwTHO=cwTHO-repmat(mean(cwTHO),CAR*par.fsDown/SLIDE,1);
    cwABD=cwABD-repmat(mean(cwABD),CAR*par.fsDown/SLIDE,1);
    for ii=1:numTap
        covTA(ii,ss)=dot(cwTHO(:,ii),cwABD(:,ii))/norm(cwTHO(:,ii))/norm(cwABD(:,ii));
    end
    
    %%State Machine
    if state(ss-1)==NOR
        rtoTHO(:,ss)=(quantile(abs((curTHO)),par.qt)./quantile(abs((preTHO)),par.qt))';%mean(abs(curTHO))/mean(abs(stdTHO));
        rtoABD(:,ss)=(quantile(abs((curABD)),par.qt)./quantile(abs((preABD)),par.qt))';%mean(abs(curABD))/mean(abs(stdABD));       
        %% Classification
        if(ss==PAR+1) %%Initialization
            Sample=[rtoABD(:,ss) rtoTHO(:,ss) rfABD(:,ss) rfTHO(:,ss)];
            fNSA(:,1) = sign(kfunNSA(svNSA,Sample,kfunargsNSA{:})'*alphaHatNSA(:) + biasNSA);
            Sample=[rtoABD(:,ss) rtoTHO(:,ss) rfABD(:,ss) rfTHO(:,ss)];
            fOC(:,1)  = sign(kfunOC(svOC,Sample,kfunargsOC{:})'*alphaHatOC(:) + biasOC);
        else
            Sample=[rtoABD(numTap,ss) rtoTHO(numTap,ss) rfABD(numTap,ss) rfTHO(numTap,ss)];
            fNSA(1:numTap-1,1)=fNSA(2:numTap,1);
            fNSA(numTap,1) = sign(kfunNSA(svNSA,Sample,kfunargsNSA{:})'*alphaHatNSA(:) + biasNSA);
            Sample=[rtoABD(numTap,ss) rtoTHO(numTap,ss) rfABD(numTap,ss) rfTHO(numTap,ss)];
            fOC(1:numTap-1,1)=fOC(2:numTap,1);
            fOC(numTap,1) = sign(kfunOC(svOC,Sample,kfunargsOC{:})'*alphaHatOC(:) + biasOC);
        end
        if(sum(fNSA)==numTap)%(sum(fNSA)>=0)%(isev==numTap)
            flagOSA=0;
            nOSA=0;
            for ff=1:numTap
                if( (covTA(ff,ss)<0)&&(fOC(ff,1)<0 ))
                   flagOSA=1;
                   nOSA=nOSA+1;
                elseif(fOC(ff,1)<0)
                    nOSA=nOSA+1;
                end
            end
            if(flagOSA==1)
                state(ss)=OSA;
                AlgoEvNum.osa=AlgoEvNum.osa+1;
            elseif(sum(fOC)==numTap)%(sum(fOC)>=0)%(vCSA>vOSA)
                state(ss)=CSA;
                AlgoEvNum.csa=AlgoEvNum.csa+1;
            elseif( (nOSA==1)&&(flagOSA==0) )
                state(ss)=CSA;
                AlgoEvNum.csa=AlgoEvNum.csa+1; 
            else
                state(ss)=OSA;
                AlgoEvNum.osa=AlgoEvNum.osa+1;
            end
            saF=1;
            tF=tF+1;
        else
            state(ss)=NOR;
            saF=0;
            tF=0;
        end       
    else
        rtoTHO(:,ss)=(quantile(abs((curTHO)),par.qt)./quantile(abs((stdTHO)),par.qt))';
        rtoABD(:,ss)=(quantile(abs((curABD)),par.qt)./quantile(abs((stdABD)),par.qt))';
        
        Sample=[rtoABD(numTap,ss) rtoTHO(numTap,ss) rfABD(numTap,ss) rfTHO(numTap,ss)];
        fNSA(1:numTap-1,1)=fNSA(2:numTap,1);
        fNSA(numTap,1) = sign(kfunNSA(svNSA,Sample,kfunargsNSA{:})'*alphaHatNSA(:) + biasNSA);
        %fNSA
        Sample=[rtoABD(numTap,ss) rtoTHO(numTap,ss) rfABD(numTap,ss) rfTHO(numTap,ss)];
        fOC(1:numTap-1,1)=fOC(2:numTap,1);
        fOC(numTap,1) = sign(kfunOC(svOC,Sample,kfunargsOC{:})'*alphaHatOC(:) + biasOC);
        if( (sum(fNSA)==(-numTap)) || (tF==120) )
            state(ss)=NOR;
            saF=0;
            tF=0;
        else
            state(ss)=state(ss-1);
            saF=0;
            tF=tF+1;
        end
    end   
end
STAGE=[STAGE';STAGE'];
STAGE=STAGE(:);
% [h]=PlotALL_v2(THO,ABD,CFlow,STAGE,rtoTHO(1,:),rtoABD(1,:),state,statePSG,par.fs,LEN,SLIDE,OSA,CSA,MSA,HYP,NOR,par.name,nev,AlgoEvNum);
feat.rt.t=rtoTHO(1,:);
feat.rt.a=rtoABD(1,:);
feat.rf.t=rfTHO(1,:);
feat.rf.a=rfABD(1,:);
feat.cov=covTA(1,:);
disp('State Machine Done!');
end
