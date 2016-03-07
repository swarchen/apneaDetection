close all 
clear all
clc
[par]=GetPars;
pnum = csvread('pnum47.csv',0,0);
len=length(pnum);
% load csa.mat

for numTap=12:12
    tic
    numTap
for pn=1:len
    close all
    par.name=num2str(pnum(pn));
    P_NUM=par.name%'1544';
%     par.groute=['./Graph/' par.name '/'];
%     mkdir(par.groute);
%     delete(['./Graph/' (P_NUM) '.pdf']);
%% SSTorNot
doSST=par.doSST;%0;
noWK=par.noWK;%1;
OSA=par.st.osa;
CSA=par.st.csa;%-1;
NOR=par.st.nor;%0;
HYP=par.st.hyp;%2;
MSA=par.st.msa;%3;

%% Read Event List
[nev,t,du]=GetEvent(P_NUM);
% load(['D:\MATLAB\Breath_Project\GrapPsgInfo\' (P_NUM) '_nev.mat']);
% load(['/home/yylin/BreathAlgo/GrapPsgInfo/' (P_NUM) '_nev.mat']);
% clearvars -except t_lightoff P_NUM doSST HI THsel FQsel noWK OSA CSA NOR HYP MSA t du nev par pnum pn Sen Acc csa osa msa hyp nor SS RR

%% Read PSG
type='all';
[ABD,THO,CFlow,STAGE]=ReadPsgFiles(P_NUM,type,par.fs);
load(['/home/yylin/BreathAlgo/RcSig1/sstTHOABD_' (P_NUM) '.mat']);
%% Parameters
fs=par.fs;%100;
fsDown=par.fsDown;%10;
downRate=par.downRate;%fs/fsDown;
t_axis=(1:length(CFlow))/fs;
T_LEN=fix(length(t_axis)/fs);
SLIDE=par.slide;%2;%1/SLIDE(sec)
CFlow=CFlow(1:downRate:end);
THO=THO(1:downRate:end);
ABD=ABD(1:downRate:end);
rcTHO=real(rcTHO)';
rcABD=real(rcABD)';
% clearvars rcTHO rcABD

%% CALL EVENT FUNCTION
[state]=GetPsgState(T_LEN,t,du,nev,par);
% load(['D:\MATLAB\Breath_Project\GrapPsgInfo\' (P_NUM) '_PsgStateInfo.mat']);
% load(['/home/yylin/BreathAlgo/GrapPsgInfo/' (P_NUM) '_PsgStateInfo.mat']);
clearvars ndata text alldata
%% Ratio Extraction
    sp=par.par+1;
    ep=min([state.len-par.fr+1,state.len-par.cov+1,state.len-par.car+1]);
    [mfTHO,mfABD,rfTHO,rfABD,fq]=freqExt3F(par,state.len,STAGE,state.psg,THO,ABD,sp,ep);    
    [covTA,cov]=covExt(par,state.len,STAGE,state.psg,THO,ABD,sp,ep,CFlow);
    [rtoTHO,rtoABD,rt]=rtExt(par,state.len,STAGE,state.psg,rcTHO,rcABD,sp,ep);
%% Ratio Analysis
    binranges =  0:0.1:3;
    islog=0;
    [rt]=rtAnalysis(rtoTHO,rtoABD,par,rt,state.psg,nev,binranges,islog);
    binranges =  0:0.01:1;
    [rf]=fqrtAnalysis(rfTHO,rfABD,par,fq,state.psg,nev,binranges);
    binranges =  -1:0.1:1;
    [cov]=covAnalysis(covTA,par,cov,state.psg,nev,binranges);
    %% If you have pooled database
%     load(['/home/yylin/BreathAlgo/feat_sstCov/' P_NUM '_covrtrf.mat']);
%     load('/home/yylin/BreathAlgo/feat_sstCov/csa_rfrtcov47.mat');
%     rt.t.csa=csa.rt.t;
%     rt.a.csa=csa.rt.a;
%     rf.t.csa=csa.rf.t;
%     rf.a.csa=csa.rf.a;
%     cov.c.csa=csa.cov;
    
    %%
    if ~isempty(rf.t.csa)
        rf.t.csa=log10(rf.t.csa);
        rf.a.csa=log10(rf.a.csa);
    end
    if ~isempty(rf.t.osa)
        rf.t.osa=log10(rf.t.osa);
        rf.a.osa=log10(rf.a.osa);
    end
    if ~isempty(rf.t.msa)
        rf.t.msa=log10(rf.t.msa);
        rf.a.msa=log10(rf.a.msa);
    end
    if ~isempty(rf.t.nor)
        rf.t.nor=log10(rf.t.nor);
        rf.a.nor=log10(rf.a.nor);
    end

%% DO SVM
trianSize=1000;
[NSAidx,svmStructNSA]=doSVM_Ar_NSA_artiGen_v4(rt,rf,par,trianSize);
[GroupOC,OCidx,svmStructOC,conf]=doSVM_ArFr_OC_artiGen_v4(rt,rf,par,trianSize,svmStructNSA);
[stateAlg,AlgoEvNum,feat]=OvNightAnl_v14(par,doSST,state.len,rcTHO,rcABD,THO,ABD,CFlow,STAGE,state.psg,nev,svmStructNSA,svmStructOC,numTap);
% save([(P_NUM) '_Tap_' num2str(numTap) '_tightCSA_all.mat'],'stateAlg','AlgoEvNum','NSAidx','OCidx','feat');
end
toc
end
