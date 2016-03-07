function [par]=GetPars

par.name='1416';
%% SSTorNot
par.doSST=0;
% par.THsel=1;
% HI=0.5:0.5:6;
par.noWK=1;
par.st.osa=1;
par.st.csa=-1;
par.st.nor=0;
par.st.hyp=2;
par.st.msa=3;
par.fs=100;
par.fsDown=4;
par.downRate=fix(par.fs/par.fsDown);
% par.t_axis=(1:length(CFlow_in))/fs;
% par.T_LEN=fix(length(t_axis)/fs);
par.slide=2;%1/SLIDE(sec)
par.cur=10;
par.pre=60;
par.fwin=10;
par.cwin=10;
par.car=par.cur*par.slide;%/SLIDE=sec%current average range
par.par=par.pre*par.slide;%/SLIDE=sec%previous average range
par.fr=par.fwin*par.slide;
par.cov=par.cwin*par.slide;
par.qt=0.95;
% par.spec.oneg=8;
% par.spec.ort=8;
% par.spec.cfq=7;
% par.spec.crt=7;
%%Hyperplane
% par.p.a=25.2152;
% par.p.b=-20.4043;
end
