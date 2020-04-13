function [t,V,I,T,Qc_cum,Qd_cum] = importingPEC_ageing_volt_concat(fname)
% fname = 'Test1126.csv';

R1 = 22-1;
% C1 = 8-1;
C1 = 8-1-1;

M = csvread(fname,R1,C1);

t = M(:,2);
V = M(:,3);
I = M(:,4);

Qc_cum = M(:,5);
Qd_cum = M(:,6);

T = M(:,7);


