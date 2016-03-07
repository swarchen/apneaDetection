function [ABD,THO,CFlow,STAGE,spo2]=ReadPsgFiles(P_NUM,type,fs)
% ABD_in= csvread(['D:\MATLAB\Breath_Project\psgdata\p' P_NUM '\p' P_NUM '_' type '\ABD.csv'],1,0);
% THO_in= csvread(['D:\MATLAB\Breath_Project\psgdata\p' P_NUM '\p' P_NUM '_' type '\THO.csv'],1,0);
% CFlow_in= csvread(['D:\MATLAB\Breath_Project\psgdata\p' P_NUM '\p' P_NUM '_' type '\CFlow.csv'],1,0);
% STAGE= csvread(['D:\MATLAB\Breath_Project\psgdata\p' P_NUM '\p' P_NUM '_' type '\STAGE.csv'],1,0);
% spo2= csvread(['D:\MATLAB\Breath_Project\psgdata\p' P_NUM '\p' P_NUM '_' type '\SpO2.csv'],1,0);
ABD_in= csvread(['/home/yylin/BreathAlgo/psgdata/p' P_NUM '/p' P_NUM '_' type '/ABD.csv'],1,0);
THO_in= csvread(['/home/yylin/BreathAlgo/psgdata/p' P_NUM '/p' P_NUM '_' type '/THO.csv'],1,0);
CFlow_in= csvread(['/home/yylin/BreathAlgo/psgdata/p' P_NUM '/p' P_NUM '_' type '/CFlow.csv'],1,0);
STAGE= csvread(['/home/yylin/BreathAlgo/psgdata/p' P_NUM '/p' P_NUM '_' type '/STAGE.csv'],1,0);
spo2= csvread(['/home/yylin/BreathAlgo/psgdata/p' P_NUM '/p' P_NUM '_' type '/SpO2.csv'],1,0);
spo2=spo2/100;

ABD = ABD_in'-2^15;
THO = THO_in'-2^15;
CFlow=CFlow_in'-2^15;
end