function [state]=GetPsgState(T_LEN,t,du,nev,par)
nOSA=nev.osa;
nCSA=nev.csa;
nMSA=nev.msa;
nHYP=nev.hyp;
t_osa=t.osa;
t_csa=t.csa;
t_msa=t.msa;
t_hyp=t.hyp;
Do=du.osa;
Dc=du.csa;
Dm=du.msa;
Dh=du.hyp;
OSA=par.st.osa;
CSA=par.st.csa;
MSA=par.st.msa;
HYP=par.st.hyp;
NOR=par.st.nor;
P_NUM=par.name;
SLIDE=par.slide;
statePSG=zeros(T_LEN*SLIDE,1);
stateHYP=zeros(T_LEN*SLIDE,1);
stateOSA=zeros(T_LEN*SLIDE,1);
stateCSA=zeros(T_LEN*SLIDE,1);
stateMSA=zeros(T_LEN*SLIDE,1);
stateNOR=ones(T_LEN*SLIDE,1);
if nOSA>0
for ii=1:length(t_osa)
    for jj=t_osa(ii)*SLIDE:(t_osa(ii)+Do(ii))*SLIDE
        statePSG(fix(jj),1)=OSA;
        stateOSA(fix(jj),1)=1;
        stateNOR(fix(jj),1)=0;
    end
end
end
if nCSA>0
for ii=1:length(t_csa)
    for jj=t_csa(ii)*SLIDE:(t_csa(ii)+Dc(ii))*SLIDE
        statePSG(fix(jj),1)=CSA;
        stateCSA(fix(jj),1)=1;
        stateNOR(fix(jj),1)=0;
    end
end
end
if nMSA>0
for ii=1:length(t_msa)
    for jj=t_msa(ii)*SLIDE:(t_msa(ii)+Dm(ii))*SLIDE
        statePSG(fix(jj),1)=MSA;
        stateMSA(fix(jj),1)=1;
        stateNOR(fix(jj),1)=0;
    end
end
end
if nHYP>0
for ii=1:length(t_hyp)
    for jj=t_hyp(ii)*SLIDE:(t_hyp(ii)+Dh(ii))*SLIDE
        statePSG(fix(jj),1)=HYP;
        stateHYP(fix(jj),1)=1;
        stateNOR(fix(jj),1)=0;
    end
end
end
LEN=min([length(statePSG),length(stateHYP),length(stateOSA),length(stateCSA),length(stateMSA),length(stateHYP)]);
state.psg=statePSG(1:LEN);
state.hyp=stateHYP(1:LEN);
state.osa=stateOSA(1:LEN);
state.csa=stateCSA(1:LEN);
state.msa=stateMSA(1:LEN);
state.nor=stateNOR(1:LEN);

% statePSG=statePSG/2+1;
state.ax=(1:length(statePSG))/2;
state.len=LEN;
end