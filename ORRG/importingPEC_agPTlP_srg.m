clear all
clc

exp11 = '1460';
exp12 = '1461';
exp21 = '1463';   Qcexp216 = 2725.130; Qdexp216 = 2721.090; tcexp216 = 4924.310; Icmxp216 = 1.0;  Icexp216 = 0.5;  Vmexp216 = 4.2;  Veexp216 = 4.2;  Vc0exp216 = 4.2;  Vd0exp216 = 4.2;  Tmexp216 = 25;
exp22 = '1465';   Qcexp22 = 2750.544;  Qdexp22 = 2747.947;  tcexp22 = 4567.400;  Icmxp22 = 1.0;   Icexp22 = 0.5;   Vmexp22 = 4.2;   Veexp22 = 4.2;   Vc0exp22 = 4.2;   Vd0exp22 = 4.2;   Tmexp22 = 25;
exp23 = '1466';
exp24 = '1467';   Qcexp24 = 2690.410;  Qdexp24 = 2687.437;  tcexp24 = 4374.400;  Icmxp24 = 1.0;   Icexp24 = 0.5;   Vmexp24 = 4.2;   Veexp24 = 4.2;   Vc0exp24 = 4.2;   Vd0exp24 = 4.2;   Tmexp24 = 25;
exp25 = '1468';   Qcexp251 = 2616.238; Qdexp251 = 2726.699; tcexp251 = 4280.210; Icmxp251 = 1.0;   Icexp251 = 0.5;  Vmexp251 = 4.2;  Veexp251 = 4.2;  Vc0exp251 = 4.2;  Vd0exp251 = 4.2;  Tmexp251 = 25;
exp26 = '1469';
exp31 = '1471';
exp32 = '1482';
exp41 = '1484';
exp42 = '1495';
exp51 = '1497';
exp52 = '1498';
exp61 = '1512';
exp62 = '1514';
exp71 = '1533';   Qcexp71 = 2476.608;                       tcexp71 = 3924.880;  Icmxp71 = 1.0;   Icexp71 = 0.5;   Vmexp71 = 4.2;   Veexp71 = 4.2;   Vc0exp71 = 4.2;                    Tmexp71 = 25;
exp72 = '1534';                        Qdexp72 = 2471.278;                                                                                                              Vd0exp72 = 4.2;
exp73 = '1535';
exp74 = '1536';
exp81 = '1541';
exp82 = '1544';
exp91 = '1547';
exp92 = '1549';   Qcexp92 = 2220.798; Qdexp92 = 2218.352; tcexp92 = 2767.183; Icmxp92 = 1.0;  Icexp92 = 0.5;  Vmexp92 = 4.4;  Veexp92 = 4.2;  Vc0exp92 = 4.2;  Vd0exp92 = 4.2;  Tmexp92 = 25;
       %'1550';
exp93 = '1552';   Qcexp93 = 2211.454; Qdexp93 = 2214.535; tcexp93 = 2684.814; Icmxp93 = 1.0;  Icexp93 = 0.5;  Vmexp93 = 4.4;  Veexp93 = 4.2;  Vc0exp93 = 4.2;  Vd0exp93 = 4.2;  Tmexp93 = 25;
       %'1553'
exp94 = '1555';   Qcexp94 = 2190.650; Qdexp94 = 2172.345; tcexp94 = 2922.230; Icmxp94 = 1.0;  Icexp94 = 0.5;  Vmexp94 = 4.4;  Veexp94 = 4.2;  Vc0exp94 = 4.2;  Vd0exp94 = 4.2;  Tmexp94 = 25;
       %'1557'
exp95 = '1560';   Qcexp95 = 2386.444; Qdexp95 = 2166.812; tcexp95 = 2975.430; Icmxp95 = 1.0;  Icexp95 = 0.5;  Vmexp95 = 4.4;  Veexp95 = 4.2;  Vc0exp95 = 4.2;  Vd0exp95 = 4.2;  Tmexp95 = 25;
       %'1570'
exp96 = '1572';   Qcexp96 = 2359.262; Qdexp96 = 2032.283; tcexp96 = 2869.279; Icmxp96 = 1.0;  Icexp96 = 0.5;  Vmexp96 = 4.4;  Veexp96 = 4.2;  Vc0exp96 = 4.2;  Vd0exp96 = 4.2;  Tmexp96 = 25;
       %'1573'
exp97 = '1575';   Qcexp97 = 2348.883; Qdexp97 = 2147.829; tcexp97 = 2924.531; Icmxp97 = 1.0;  Icexp97 = 0.5;  Vmexp97 = 4.4;  Veexp97 = 4.2;  Vc0exp97 = 4.2;  Vd0exp97 = 4.2;  Tmexp97 = 25;
       %'1576'
exp98 = '1578';   Qcexp98 = 2272.674; Qdexp98 = 2172.133; tcexp98 = 2823.977; Icmxp98 = 1.0;  Icexp98 = 0.5;  Vmexp98 = 4.4;  Veexp98 = 4.2;  Vc0exp98 = 4.2;  Vd0exp98 = 4.2;  Tmexp98 = 25;
       %'1579'
exp99 = '1581';   Qcexp99 = 2282.239; Qdexp99 = 1992.952; tcexp99 = 2689.955; Icmxp99 = 1.0;  Icexp99 = 0.5;  Vmexp99 = 4.4;  Veexp99 = 4.2;  Vc0exp99 = 4.2;  Vd0exp99 = 4.2;  Tmexp99 = 25;
       %'1582'
exp01 = '1583';

expM = {exp11; [exp21; exp22; exp23; exp24; exp25]; exp31; exp41; exp51; exp61; [exp71; exp73]; exp81; [exp91; exp94]};
expP = {exp12; exp26; exp32; exp42; exp52; exp62; exp74; exp82; exp01};% exp4; exp6; exp8; exp10];
expN = {11; [5 1 1 1 2]; 11; 11; 11; 11; [3 7]; 11; [2 1]};% 11; 11; 11; 11];
% spec = {0; [6 1 1 1 1]};
expR = {0; [Qcexp216; Qdexp216; tcexp216; Icmxp216;  Icexp216;  Vmexp216;  Veexp216;  Vc0exp216;  Vd0exp216;  Tmexp216]; ...
    [Qcexp22;  Qdexp22;  tcexp22;  Icmxp22;  Icexp22;   Vmexp22;   Veexp22;   Vc0exp22;   Vd0exp22;   Tmexp22]; ...
    0; [Qcexp24;  Qdexp24;  tcexp24;  Icmxp24;  Icexp24;   Vmexp24;   Veexp24;   Vc0exp24;   Vd0exp24;   Tmexp24]; ...
    [Qcexp251; Qdexp251; tcexp251; Icmxp251; Icexp251;  Vmexp251;  Veexp251;  Vc0exp251;  Vd0exp251;  Tmexp251]};
n_exp = 9;

ncyc_exp = 11 + 1;
vcyc = 1:2*ncyc_exp*n_exp;
Mcyc = reshape(vcyc,2*ncyc_exp,n_exp)';
% break

% count = 26;
% count = str2double(exp1(end-1:end));
% Qc_a_tot = zeros(1,[];
% Qd_a_tot = [];
neven = 2:2:2*ncyc_exp-1;
nodd = 1:2:2*ncyc_exp-2;
% L = strsplit(sprintf('%c\n','1':'5'));  % Letter Labels
fg0 = 80;
% break

%%
j=1;
expM0 = expM{j};
expN0 = expN{j};
expP0 = expP{j};
% spec0 = spec{j};
count = str2double(expM0(1,:));
fname_a = sprintf('Test%d.csv',count);
[Qc_a,Qd_a, Vdc_a,Vdd_a, tc_a,td_a, Icm_a,Ic_a,Id_a,Vcm_a,Vce_a, Vc0_a,Vd0_a, Tcm_a] = importingPEC_ageing(fname_a, expN0(1));
Pc_a = Vdc_a * 5.87;
Pd_a = Vdd_a * 5.87/2;

vidx = 1:length(Qc_a);
figure(fg0)
plot(Mcyc(j,nodd(vidx)), Qc_a/1e3,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;
plot(Mcyc(j,neven(vidx)),Qd_a/1e3,'og','MarkerFaceColor','g','markersize',6);hold on;axis tight;%title('Capacity fade')
figure(fg0+2)
plot(Mcyc(j,nodd(2:vidx(end))),tc_a(2:end)/3600,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
figure(fg0+3)
plot(Mcyc(j,nodd(2:vidx(end))),Icm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Average I')
plot(Mcyc(j,nodd(2:vidx(end))),Ic_a(2:end),'om','markersize',6);hold on;axis tight;%title('Average I')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+4)
plot(Mcyc(j,nodd(2:vidx(end))),Vcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(2:vidx(end))),Vce_a(2:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+5)
plot(Mcyc(j,nodd(2:vidx(end))),Tcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max T')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[20 70],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+6)
plot(Mcyc(j,nodd(1:vidx(end))),Vc0_a(1:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(1:vidx(end))),Vd0_a(1:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+7)
plot(Mcyc(j,nodd(2:vidx(end))),(Qc_a(2:end)/1e3)./(tc_a(2:end)/3600),'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
vV0 = [Vc0_a(1:end); Vd0_a(1:end)];
% break

count = str2double(expP0(1,:));
fname_p = sprintf('Test%d.csv',count);
[Qc_p,Qd_p, Vdc_p,Vdd_p, tc_p,td_p] = importingPEC_PT(fname_p);
Pc_p = Vdc_p * 5.87/2;
Pd_p = Vdd_p * 5.87/2;

figure(fg0)
plot(Mcyc(j,end-1),Qc_p   /1e3,'or','linewidth',2,'markersize',12);hold on;axis tight;
plot(Mcyc(j,end),  Qd_p(2)/1e3,'og','linewidth',2,'markersize',12);hold on;axis tight;
% plot((Mcyc(j,end)+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
% break
%%
j=2;    k=1;
expM0 = expM{j};
expN0 = expN{j};
expP0 = expP{j};

count = str2double(expM0(1,:));
fname_a = sprintf('Test%d.csv',count);
[Qc_a,Qd_a, Vdc_a,Vdd_a, tc_a,td_a, Icm_a,Ic_a,Id_a,Vcm_a,Vce_a, Vc0_a,Vd0_a, Tcm_a] = importingPEC_ageing(fname_a, expN0(1));
Pc_a = Vdc_a * 5.87;
Pd_a = Vdd_a * 5.87/2;

% k=k+1;
Qc_a = [Qc_a Qcexp216]; Qd_a = [Qd_a Qdexp216]; Vdc_a = [Vdc_a ]; Vdd_a = [Vdd_a ]; tc_a = [tc_a tcexp216];
td_a = [td_a ]; Icm_a = [Icm_a Icmxp216]; Ic_a = [Ic_a Icexp216]; Id_a = [Id_a ];   Vcm_a = [Vcm_a Vmexp216]; 
Vce_a = [Vce_a Veexp216]; Vc0_a = [Vc0_a Vc0exp216]; Vd0_a = [Vd0_a Vd0exp216]; Tcm_a = [Tcm_a Tmexp216];

k=k+1;
Qc_a = [Qc_a Qcexp22]; Qd_a = [Qd_a Qdexp22]; Vdc_a = [Vdc_a ]; Vdd_a = [Vdd_a ]; tc_a = [tc_a tcexp22];
td_a = [td_a ]; Icm_a = [Icm_a Icmxp22]; Ic_a = [Ic_a Icexp22]; Id_a = [Id_a ];   Vcm_a = [Vcm_a Vmexp22]; 
Vce_a = [Vce_a Veexp22]; Vc0_a = [Vc0_a Vc0exp22]; Vd0_a = [Vd0_a Vd0exp22]; Tcm_a = [Tcm_a Tmexp22];

k=k+1;
count = str2double(expM0(k,:));
fname_a = sprintf('Test%d.csv',count);
[Qc_a2,Qd_a2, Vdc_a2,Vdd_a2, tc_a2,td_a2, Icm_a2,Ic_a2,Id_a2,Vcm_a2,Vce_a2, Vc0_a2,Vd0_a2, Tcm_a2] = importingPEC_ageing(fname_a, expN0(k));
Qc_a = [Qc_a Qc_a2]; Qd_a = [Qd_a Qd_a2]; Vdc_a = [Vdc_a Vdc_a2]; Vdd_a = [Vdd_a Vdd_a2]; tc_a = [tc_a tc_a2];
td_a = [td_a td_a2]; Icm_a = [Icm_a Icm_a2]; Ic_a = [Ic_a Ic_a2]; Id_a = [Id_a Id_a2];   Vcm_a = [Vcm_a Vcm_a2]; 
Vce_a = [Vce_a Vce_a2]; Vc0_a = [Vc0_a Vc0_a2]; Vd0_a = [Vd0_a Vd0_a2]; Tcm_a = [Tcm_a Tcm_a2];
% k
% expN0(k)
% expM0(k,:)
% Qc_a
% pause

k=k+1;
Qc_a = [Qc_a Qcexp24]; Qd_a = [Qd_a Qdexp24]; Vdc_a = [Vdc_a ]; Vdd_a = [Vdd_a ]; tc_a = [tc_a tcexp24];
td_a = [td_a ]; Icm_a = [Icm_a Icmxp24]; Ic_a = [Ic_a Icexp24]; Id_a = [Id_a ];   Vcm_a = [Vcm_a Vmexp24]; 
Vce_a = [Vce_a Veexp24]; Vc0_a = [Vc0_a Vc0exp24]; Vd0_a = [Vd0_a Vd0exp24]; Tcm_a = [Tcm_a Tmexp24];
% Qc_a
% pause

k=k+1;
Qc_a = [Qc_a Qcexp251]; Qd_a = [Qd_a Qdexp251]; Vdc_a = [Vdc_a ]; Vdd_a = [Vdd_a ]; tc_a = [tc_a tcexp251];
td_a = [td_a ]; Icm_a = [Icm_a Icmxp251]; Ic_a = [Ic_a Icexp251]; Id_a = [Id_a ];   Vcm_a = [Vcm_a Vmexp251]; 
Vce_a = [Vce_a Veexp251]; Vc0_a = [Vc0_a Vc0exp251]; Vd0_a = [Vd0_a Vd0exp251]; Tcm_a = [Tcm_a Tmexp251];
% Qc_a
% pause

count = str2double(expM0(k,:));
fname_a = sprintf('Test%d.csv',count);
[Qc_a3,Qd_a3, Vdc_a3,Vdd_a3, tc_a3,td_a3, Icm_a3,Ic_a3,Id_a3,Vcm_a3,Vce_a3,Vc0_a3,Vd0_a3,Tcm_a3] = importingPEC_ageing2(fname_a, expN0(k), 2);
Qc_a = [Qc_a Qc_a3(2)]; Qd_a = [Qd_a Qd_a3(2)]; Vdc_a = [Vdc_a Vdc_a3(2)]; Vdd_a = [Vdd_a Vdd_a3(2)]; tc_a = [tc_a tc_a3(2)];
td_a = [td_a td_a3(2)]; Icm_a = [Icm_a Icm_a3(2)]; Ic_a = [Ic_a Ic_a3(2)]; Id_a = [Id_a Id_a3(2)];   Vcm_a = [Vcm_a Vcm_a3(2)]; 
Vce_a = [Vce_a Vce_a3(2)]; Vc0_a = [Vc0_a Vc0_a3(2)]; Vd0_a = [Vd0_a Vd0_a3(2)]; Tcm_a = [Tcm_a Tcm_a3(2)];
% k
% expN0(k)
% expM0(k,:)
% Qc_a
% pause

vidx = 1:length(Qc_a);
figure(fg0)
plot(Mcyc(j,nodd(vidx)), Qc_a/1e3,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;
plot(Mcyc(j,neven(vidx)),Qd_a/1e3,'og','MarkerFaceColor','g','markersize',6);hold on;axis tight;%title('Capacity fade')
figure(fg0+2)
plot(Mcyc(j,nodd(2:vidx(end))),tc_a(2:end)/3600,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
figure(fg0+3)
plot(Mcyc(j,nodd(2:vidx(end))),Icm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Average I')
plot(Mcyc(j,nodd(2:vidx(end))),Ic_a(2:end),'om','markersize',6);hold on;axis tight;%title('Average I')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+4)
plot(Mcyc(j,nodd(2:vidx(end))),Vcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(2:vidx(end))),Vce_a(2:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+5)
plot(Mcyc(j,nodd(2:vidx(end))),Tcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max T')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[20 70],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+6)
plot(Mcyc(j,nodd(1:vidx(end))),Vc0_a(1:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(1:vidx(end))),Vd0_a(1:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+7)
plot(Mcyc(j,nodd(2:vidx(end))),(Qc_a(2:end)/1e3)./(tc_a(2:end)/3600),'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
vV0 = [vV0 [Vc0_a(1:end); Vd0_a(1:end)]];

count = str2double(expP0(1,:));
fname_p = sprintf('Test%d.csv',count);
[Qc_p,Qd_p, Vdc_p,Vdd_p, tc_p,td_p] = importingPEC_PT(fname_p);
Pc_p = Vdc_p * 5.87/2;
Pd_p = Vdd_p * 5.87/2;

figure(fg0)
plot(Mcyc(j,end-1),Qc_p   /1e3,'or','linewidth',2,'markersize',12);hold on;axis tight;
plot(Mcyc(j,end),  Qd_p(2)/1e3,'og','linewidth',2,'markersize',12);hold on;axis tight;
% plot((Mcyc(j,end)+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
% break
%%
j=3;    k=1;
expM0 = expM{j};
expN0 = expN{j};
expP0 = expP{j};
% % % spec0 = spec{j};
count = str2double(expM0(k,:));
fname_a = sprintf('Test%d.csv',count);
[Qc_a,Qd_a, Vdc_a,Vdd_a, tc_a,td_a, Icm_a,Ic_a,Id_a,Vcm_a,Vce_a, Vc0_a,Vd0_a, Tcm_a] = importingPEC_ageing(fname_a, expN0(k));
Pc_a = Vdc_a * 5.87;
Pd_a = Vdd_a * 5.87/2;

vidx = 1:length(Qc_a);
figure(fg0)
plot(Mcyc(j,nodd(vidx)), Qc_a/1e3,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;
plot(Mcyc(j,neven(vidx)),Qd_a/1e3,'og','MarkerFaceColor','g','markersize',6);hold on;axis tight;%title('Capacity fade')
figure(fg0+2)
plot(Mcyc(j,nodd(2:vidx(end))),tc_a(2:end)/3600,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
figure(fg0+3)
plot(Mcyc(j,nodd(2:vidx(end))),Icm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Average I')
plot(Mcyc(j,nodd(2:vidx(end))),Ic_a(2:end),'om','markersize',6);hold on;axis tight;%title('Average I')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+4)
plot(Mcyc(j,nodd(2:vidx(end))),Vcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(2:vidx(end))),Vce_a(2:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+5)
plot(Mcyc(j,nodd(2:vidx(end))),Tcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max T')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[20 70],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+6)
plot(Mcyc(j,nodd(1:vidx(end))),Vc0_a(1:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(1:vidx(end))),Vd0_a(1:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
% break
figure(fg0+7)
plot(Mcyc(j,nodd(2:vidx(end))),(Qc_a(2:end)/1e3)./(tc_a(2:end)/3600),'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
vV0 = [vV0 [Vc0_a(1:end); Vd0_a(1:end)]];

count = str2double(expP0(k,:));
fname_p = sprintf('Test%d.csv',count);
[Qc_p,Qd_p, Vdc_p,Vdd_p, tc_p,td_p] = importingPEC_PT(fname_p);
Pc_p = Vdc_p * 5.87/2;
Pd_p = Vdd_p * 5.87/2;

figure(fg0)
plot(Mcyc(j,end-1),Qc_p   /1e3,'or','linewidth',2,'markersize',12);hold on;axis tight;
plot(Mcyc(j,end),  Qd_p(2)/1e3,'og','linewidth',2,'markersize',12);hold on;axis tight;
% plot((Mcyc(j,end)+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
%%
j=4;    k=1;
expM0 = expM{j};
expN0 = expN{j};
expP0 = expP{j};
% % % spec0 = spec{j};
count = str2double(expM0(k,:));
fname_a = sprintf('Test%d.csv',count);
[Qc_a,Qd_a, Vdc_a,Vdd_a, tc_a,td_a, Icm_a,Ic_a,Id_a,Vcm_a,Vce_a, Vc0_a,Vd0_a, Tcm_a] = importingPEC_ageing(fname_a, expN0(k));
Pc_a = Vdc_a * 5.87;
Pd_a = Vdd_a * 5.87/2;

vidx = 1:length(Qc_a);
figure(fg0)
plot(Mcyc(j,nodd(vidx)), Qc_a/1e3,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;
plot(Mcyc(j,neven(vidx)),Qd_a/1e3,'og','MarkerFaceColor','g','markersize',6);hold on;axis tight;%title('Capacity fade')
figure(fg0+2)
plot(Mcyc(j,nodd(2:vidx(end))),tc_a(2:end)/3600,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
figure(fg0+3)
plot(Mcyc(j,nodd(2:vidx(end))),Icm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Average I')
plot(Mcyc(j,nodd(2:vidx(end))),Ic_a(2:end),'om','markersize',6);hold on;axis tight;%title('Average I')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+4)
plot(Mcyc(j,nodd(2:vidx(end))),Vcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(2:vidx(end))),Vce_a(2:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+5)
plot(Mcyc(j,nodd(2:vidx(end))),Tcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max T')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[20 70],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+6)
plot(Mcyc(j,nodd(1:vidx(end))),Vc0_a(1:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(1:vidx(end))),Vd0_a(1:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
% break
figure(fg0+7)
plot(Mcyc(j,nodd(2:vidx(end))),(Qc_a(2:end)/1e3)./(tc_a(2:end)/3600),'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
vV0 = [vV0 [Vc0_a(1:end); Vd0_a(1:end)]];

count = str2double(expP0(k,:));
fname_p = sprintf('Test%d.csv',count);
[Qc_p,Qd_p, Vdc_p,Vdd_p, tc_p,td_p] = importingPEC_PT(fname_p);
Pc_p = Vdc_p * 5.87/2;
Pd_p = Vdd_p * 5.87/2;

figure(fg0)
plot(Mcyc(j,end-1),Qc_p   /1e3,'or','linewidth',2,'markersize',12);hold on;axis tight;
plot(Mcyc(j,end),  Qd_p(2)/1e3,'og','linewidth',2,'markersize',12);hold on;axis tight;
% plot((Mcyc(j,end)+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
%%
j=5;    k=1;
expM0 = expM{j};
expN0 = expN{j};
expP0 = expP{j};
% % % spec0 = spec{j};
count = str2double(expM0(k,:));
fname_a = sprintf('Test%d.csv',count);
[Qc_a,Qd_a, Vdc_a,Vdd_a, tc_a,td_a, Icm_a,Ic_a,Id_a,Vcm_a,Vce_a, Vc0_a,Vd0_a, Tcm_a] = importingPEC_ageing(fname_a, expN0(k));
Pc_a = Vdc_a * 5.87;
Pd_a = Vdd_a * 5.87/2;

vidx = 1:length(Qc_a);
figure(fg0)
plot(Mcyc(j,nodd(vidx)), Qc_a/1e3,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;
plot(Mcyc(j,neven(vidx)),Qd_a/1e3,'og','MarkerFaceColor','g','markersize',6);hold on;axis tight;%title('Capacity fade')
figure(fg0+2)
plot(Mcyc(j,nodd(2:vidx(end))),tc_a(2:end)/3600,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
figure(fg0+3)
plot(Mcyc(j,nodd(2:vidx(end))),Icm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Average I')
plot(Mcyc(j,nodd(2:vidx(end))),Ic_a(2:end),'om','markersize',6);hold on;axis tight;%title('Average I')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+4)
plot(Mcyc(j,nodd(2:vidx(end))),Vcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(2:vidx(end))),Vce_a(2:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+5)
plot(Mcyc(j,nodd(2:vidx(end))),Tcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max T')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[20 70],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+6)
plot(Mcyc(j,nodd(1:vidx(end))),Vc0_a(1:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(1:vidx(end))),Vd0_a(1:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
% break
figure(fg0+7)
plot(Mcyc(j,nodd(2:vidx(end))),(Qc_a(2:end)/1e3)./(tc_a(2:end)/3600),'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
vV0 = [vV0 [Vc0_a(1:end); Vd0_a(1:end)]];

count = str2double(expP0(k,:));
fname_p = sprintf('Test%d.csv',count);
[Qc_p,Qd_p, Vdc_p,Vdd_p, tc_p,td_p] = importingPEC_PT(fname_p);
Pc_p = Vdc_p * 5.87/2;
Pd_p = Vdd_p * 5.87/2;

figure(fg0)
plot(Mcyc(j,end-1),Qc_p   /1e3,'or','linewidth',2,'markersize',12);hold on;axis tight;
plot(Mcyc(j,end),  Qd_p(2)/1e3,'og','linewidth',2,'markersize',12);hold on;axis tight;
% plot((Mcyc(j,end)+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
% break
%% Long Performance Test
fname_l = 'Test1507.csv';
[Qc_l,Qd_l, Vdc_l,Vdd_l, tc_l,td_l] = importingPEC_longPbb(fname_l);

figure(fg0)
plot(Mcyc(j,end)+[1 3],Qc_l/1e3,'pr',Mcyc(j,end)+[2 4],Qd_l/1e3,'pg','markersize',24);hold on;axis tight;%title('Capacity fade')
%%
j=6;    k=1;
expM0 = expM{j};
expN0 = expN{j};
expP0 = expP{j};
% % % spec0 = spec{j};
count = str2double(expM0(k,:));
fname_a = sprintf('Test%d.csv',count);
[Qc_a,Qd_a, Vdc_a,Vdd_a, tc_a,td_a, Icm_a,Ic_a,Id_a,Vcm_a,Vce_a, Vc0_a,Vd0_a, Tcm_a] = importingPEC_ageing(fname_a, expN0(k));
Pc_a = Vdc_a * 5.87;
Pd_a = Vdd_a * 5.87/2;

vidx = 1:length(Qc_a);
figure(fg0)
plot(Mcyc(j,nodd(vidx)), Qc_a/1e3,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight
plot(Mcyc(j,neven(vidx)),Qd_a/1e3,'og','MarkerFaceColor','g','markersize',6);hold on;axis tight;%title('Capacity fade')
figure(fg0+2)
plot(Mcyc(j,nodd(2:vidx(end))),tc_a(2:end)/3600,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
figure(fg0+3)
plot(Mcyc(j,nodd(2:vidx(end))),Icm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Average I')
plot(Mcyc(j,nodd(2:vidx(end))),Ic_a(2:end),'om','markersize',6);hold on;axis tight;%title('Average I')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+4)
plot(Mcyc(j,nodd(2:vidx(end))),Vcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(2:vidx(end))),Vce_a(2:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+5)
plot(Mcyc(j,nodd(2:vidx(end))),Tcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max T')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[20 70],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+6)
plot(Mcyc(j,nodd(1:vidx(end))),Vc0_a(1:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(1:vidx(end))),Vd0_a(1:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
% break
figure(fg0+7)
plot(Mcyc(j,nodd(2:vidx(end))),(Qc_a(2:end)/1e3)./(tc_a(2:end)/3600),'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')

vV0 = [vV0 [Vc0_a(1:end); Vd0_a(1:end)]];

count = str2double(expP0(k,:));
fname_p = sprintf('Test%d.csv',count);
[Qc_p,Qd_p, Vdc_p,Vdd_p, tc_p,td_p] = importingPEC_PT(fname_p);
Pc_p = Vdc_p * 5.87/2;
Pd_p = Vdd_p * 5.87/2;

figure(fg0)
plot(Mcyc(j,end-1),Qc_p/1e3   ,'or','linewidth',2,'markersize',12);hold on;axis tight;
plot(Mcyc(j,end),  Qd_p(2)/1e3,'og','linewidth',2,'markersize',12);hold on;axis tight;
% plot((Mcyc(j,end)+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
%%
j=7;    k=1;
expM0 = expM{j};
expN0 = expN{j};
expP0 = expP{j};
% % % spec0 = spec{j};
count = str2double(expM0(k,:));
fname_a = sprintf('Test%d.csv',count);
[Qc_a,Qd_a, Vdc_a,Vdd_a, tc_a,td_a, Icm_a,Ic_a,Id_a,Vcm_a,Vce_a, Vc0_a,Vd0_a, Tcm_a] = importingPEC_ageing(fname_a, expN0(k));
Pc_a = Vdc_a * 5.87;
Pd_a = Vdd_a * 5.87/2;

% k=k+1;
Qc_a = [Qc_a Qcexp71]; Qd_a = [Qd_a Qdexp72]; Vdc_a = [Vdc_a ]; Vdd_a = [Vdd_a ]; tc_a = [tc_a tcexp71];
td_a = [td_a ]; Icm_a = [Icm_a Icmxp71]; Ic_a = [Ic_a Icexp71]; Id_a = [Id_a ];   Vcm_a = [Vcm_a Vmexp71]; 
Vce_a = [Vce_a Veexp71]; Vc0_a = [Vc0_a Vc0exp71]; Vd0_a = [Vd0_a Vd0exp72]; Tcm_a = [Tcm_a Tmexp71];

k=k+1;
count = str2double(expM0(k,:));
fname_a = sprintf('Test%d.csv',count);
[Qc_a2,Qd_a2, Vdc_a2,Vdd_a2, tc_a2,td_a2, Icm_a2,Ic_a2,Id_a2,Vcm_a2,Vce_a2, Vc0_a2,Vd0_a2, Tcm_a2] = importingPEC_ageing(fname_a, expN0(k));
Qc_a = [Qc_a Qc_a2]; Qd_a = [Qd_a Qd_a2]; Vdc_a = [Vdc_a Vdc_a2]; Vdd_a = [Vdd_a Vdd_a2]; tc_a = [tc_a tc_a2];
td_a = [td_a td_a2]; Icm_a = [Icm_a Icm_a2]; Ic_a = [Ic_a Ic_a2]; Id_a = [Id_a Id_a2];   Vcm_a = [Vcm_a Vcm_a2]; 
Vce_a = [Vce_a Vce_a2]; Vc0_a = [Vc0_a Vc0_a2]; Vd0_a = [Vd0_a Vd0_a2]; Tcm_a = [Tcm_a Tcm_a2];

vidx = 1:length(Qc_a);
figure(fg0)
plot(Mcyc(j,nodd(vidx)), Qc_a/1e3,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;
plot(Mcyc(j,neven(vidx)),Qd_a/1e3,'og','MarkerFaceColor','g','markersize',6);hold on;axis tight;%title('Capacity fade')
figure(fg0+2)
plot(Mcyc(j,nodd(2:vidx(end))),tc_a(2:end)/3600,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
figure(fg0+3)
plot(Mcyc(j,nodd(2:vidx(end))),Icm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Average I')
plot(Mcyc(j,nodd(2:vidx(end))),Ic_a(2:end),'om','markersize',6);hold on;axis tight;%title('Average I')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+4)
plot(Mcyc(j,nodd(2:vidx(end))),Vcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(2:vidx(end))),Vce_a(2:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+5)
plot(Mcyc(j,nodd(2:vidx(end))),Tcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max T')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[20 70],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+6)
plot(Mcyc(j,nodd(1:vidx(end))),Vc0_a(1:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(1:vidx(end))),Vd0_a(1:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
% break
figure(fg0+7)
plot(Mcyc(j,nodd(2:vidx(end))),(Qc_a(2:end)/1e3)./(tc_a(2:end)/3600),'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
vV0 = [vV0 [Vc0_a(1:end); Vd0_a(1:end)]];

count = str2double(expP0(1,:));
fname_p = sprintf('Test%d.csv',count);
[Qc_p,Qd_p, Vdc_p,Vdd_p, tc_p,td_p] = importingPEC_PT(fname_p);
Pc_p = Vdc_p * 5.87/2;
Pd_p = Vdd_p * 5.87/2;

figure(fg0)
plot(Mcyc(j,end-1),Qc_p/1e3,   'or','linewidth',2,'markersize',12);hold on;axis tight;
plot(Mcyc(j,end),  Qd_p(2)/1e3,'og','linewidth',2,'markersize',12);hold on;axis tight;
% plot((Mcyc(j,end)+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
%%
j=8;    k=1;
expM0 = expM{j};
expN0 = expN{j};
expP0 = expP{j};
% % % spec0 = spec{j};
count = str2double(expM0(k,:));
fname_a = sprintf('Test%d.csv',count);
[Qc_a,Qd_a, Vdc_a,Vdd_a, tc_a,td_a, Icm_a,Ic_a,Id_a,Vcm_a,Vce_a, Vc0_a,Vd0_a, Tcm_a] = importingPEC_ageing(fname_a, expN0(k));
Pc_a = Vdc_a * 5.87;
Pd_a = Vdd_a * 5.87/2;

vidx = 1:length(Qc_a);
figure(fg0)
plot(Mcyc(j,nodd(vidx)), Qc_a/1e3,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;
plot(Mcyc(j,neven(vidx)),Qd_a/1e3,'og','MarkerFaceColor','g','markersize',6);hold on;axis tight;%title('Capacity fade')
figure(fg0+2)
plot(Mcyc(j,nodd(2:vidx(end))),tc_a(2:end)/3600,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
figure(fg0+3)
plot(Mcyc(j,nodd(2:vidx(end))),Icm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Average I')
plot(Mcyc(j,nodd(2:vidx(end))),Ic_a(2:end),'om','markersize',6);hold on;axis tight;%title('Average I')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+4)
plot(Mcyc(j,nodd(2:vidx(end))),Vcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(2:vidx(end))),Vce_a(2:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+5)
plot(Mcyc(j,nodd(2:vidx(end))),Tcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max T')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[20 70],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+6)
plot(Mcyc(j,nodd(1:vidx(end))),Vc0_a(1:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(1:vidx(end))),Vd0_a(1:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
% break
figure(fg0+7)
plot(Mcyc(j,nodd(2:vidx(end))),(Qc_a(2:end)/1e3)./(tc_a(2:end)/3600),'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
vV0 = [vV0 [Vc0_a(1:end); Vd0_a(1:end)]];

count = str2double(expP0(k,:));
fname_p = sprintf('Test%d.csv',count);
[Qc_p,Qd_p, Vdc_p,Vdd_p, tc_p,td_p] = importingPEC_PT(fname_p);
Pc_p = Vdc_p * 5.87/2;
Pd_p = Vdd_p * 5.87/2;

figure(fg0)
plot(Mcyc(j,end-1),Qc_p/1e3,   'or','linewidth',2,'markersize',12);hold on;axis tight;
plot(Mcyc(j,end),  Qd_p(2)/1e3,'og','linewidth',2,'markersize',12);hold on;axis tight;
% plot((Mcyc(j,end)+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
%%
j=9;    k=1;
expM0 = expM{j};
expN0 = expN{j};
expP0 = expP{j};
% % % spec0 = spec{j};
count = str2double(expM0(k,:));
fname_a = sprintf('Test%d.csv',count);
[Qc_a,Qd_a, Vdc_a,Vdd_a, tc_a,td_a, Icm_a,Ic_a,Id_a,Vcm_a,Vce_a, Vc0_a,Vd0_a, Tcm_a] = importingPEC_ageing(fname_a, expN0(k));
Pc_a = Vdc_a * 5.87;
Pd_a = Vdd_a * 5.87/2;

% k=k+1;
Qc_a = [Qc_a Qcexp92]; Qd_a = [Qd_a Qdexp92]; Vdc_a = [Vdc_a ]; Vdd_a = [Vdd_a ]; tc_a = [tc_a tcexp92];
td_a = [td_a ]; Icm_a = [Icm_a Icmxp92]; Ic_a = [Ic_a Icexp92]; Id_a = [Id_a ];   Vcm_a = [Vcm_a Vmexp92]; 
Vce_a = [Vce_a Veexp92]; Vc0_a = [Vc0_a Vc0exp92]; Vd0_a = [Vd0_a Vd0exp92]; Tcm_a = [Tcm_a Tmexp92];

% k=k+1;
Qc_a = [Qc_a Qcexp93]; Qd_a = [Qd_a Qdexp93]; Vdc_a = [Vdc_a ]; Vdd_a = [Vdd_a ]; tc_a = [tc_a tcexp93];
td_a = [td_a ]; Icm_a = [Icm_a Icmxp93]; Ic_a = [Ic_a Icexp93]; Id_a = [Id_a ];   Vcm_a = [Vcm_a Vmexp93]; 
Vce_a = [Vce_a Veexp93]; Vc0_a = [Vc0_a Vc0exp93]; Vd0_a = [Vd0_a Vd0exp93]; Tcm_a = [Tcm_a Tmexp93];

k=k+1;
count = str2double(expM0(k,:));
fname_a = sprintf('Test%d.csv',count);
[Qc_a2,Qd_a2, Vdc_a2,Vdd_a2, tc_a2,td_a2, Icm_a2,Ic_a2,Id_a2,Vcm_a2,Vce_a2, Vc0_a2,Vd0_a2, Tcm_a2] = importingPEC_ageing(fname_a, expN0(k));
Qc_a = [Qc_a Qc_a2]; Qd_a = [Qd_a Qd_a2]; Vdc_a = [Vdc_a Vdc_a2]; Vdd_a = [Vdd_a Vdd_a2]; tc_a = [tc_a tc_a2];
td_a = [td_a td_a2]; Icm_a = [Icm_a Icm_a2]; Ic_a = [Ic_a Ic_a2]; Id_a = [Id_a Id_a2];   Vcm_a = [Vcm_a Vcm_a2]; 
Vce_a = [Vce_a Vce_a2]; Vc0_a = [Vc0_a Vc0_a2]; Vd0_a = [Vd0_a Vd0_a2]; Tcm_a = [Tcm_a Tcm_a2];

% k=k+1;
Qc_a = [Qc_a Qcexp94]; Qd_a = [Qd_a Qdexp94]; Vdc_a = [Vdc_a ]; Vdd_a = [Vdd_a ]; tc_a = [tc_a tcexp94];
td_a = [td_a ]; Icm_a = [Icm_a Icmxp94]; Ic_a = [Ic_a Icexp94]; Id_a = [Id_a ];   Vcm_a = [Vcm_a Vmexp94]; 
Vce_a = [Vce_a Veexp94]; Vc0_a = [Vc0_a Vc0exp94]; Vd0_a = [Vd0_a Vd0exp94]; Tcm_a = [Tcm_a Tmexp94];

% k=k+1;
Qc_a = [Qc_a Qcexp95]; Qd_a = [Qd_a Qdexp95]; Vdc_a = [Vdc_a ]; Vdd_a = [Vdd_a ]; tc_a = [tc_a tcexp95];
td_a = [td_a ]; Icm_a = [Icm_a Icmxp95]; Ic_a = [Ic_a Icexp95]; Id_a = [Id_a ];   Vcm_a = [Vcm_a Vmexp95]; 
Vce_a = [Vce_a Veexp95]; Vc0_a = [Vc0_a Vc0exp95]; Vd0_a = [Vd0_a Vd0exp95]; Tcm_a = [Tcm_a Tmexp95];

% k=k+1;
Qc_a = [Qc_a Qcexp96]; Qd_a = [Qd_a Qdexp96]; Vdc_a = [Vdc_a ]; Vdd_a = [Vdd_a ]; tc_a = [tc_a tcexp96];
td_a = [td_a ]; Icm_a = [Icm_a Icmxp96]; Ic_a = [Ic_a Icexp96]; Id_a = [Id_a ];   Vcm_a = [Vcm_a Vmexp96]; 
Vce_a = [Vce_a Veexp96]; Vc0_a = [Vc0_a Vc0exp96]; Vd0_a = [Vd0_a Vd0exp96]; Tcm_a = [Tcm_a Tmexp96];

% k=k+1;
Qc_a = [Qc_a Qcexp97]; Qd_a = [Qd_a Qdexp97]; Vdc_a = [Vdc_a ]; Vdd_a = [Vdd_a ]; tc_a = [tc_a tcexp97];
td_a = [td_a ]; Icm_a = [Icm_a Icmxp97]; Ic_a = [Ic_a Icexp97]; Id_a = [Id_a ];   Vcm_a = [Vcm_a Vmexp97]; 
Vce_a = [Vce_a Veexp97]; Vc0_a = [Vc0_a Vc0exp97]; Vd0_a = [Vd0_a Vd0exp97]; Tcm_a = [Tcm_a Tmexp97];

% k=k+1;
Qc_a = [Qc_a Qcexp98]; Qd_a = [Qd_a Qdexp98]; Vdc_a = [Vdc_a ]; Vdd_a = [Vdd_a ]; tc_a = [tc_a tcexp98];
td_a = [td_a ]; Icm_a = [Icm_a Icmxp98]; Ic_a = [Ic_a Icexp98]; Id_a = [Id_a ];   Vcm_a = [Vcm_a Vmexp98]; 
Vce_a = [Vce_a Veexp98]; Vc0_a = [Vc0_a Vc0exp98]; Vd0_a = [Vd0_a Vd0exp98]; Tcm_a = [Tcm_a Tmexp98];

% k=k+1;
Qc_a = [Qc_a Qcexp99]; Qd_a = [Qd_a Qdexp99]; Vdc_a = [Vdc_a ]; Vdd_a = [Vdd_a ]; tc_a = [tc_a tcexp99];
td_a = [td_a ]; Icm_a = [Icm_a Icmxp99]; Ic_a = [Ic_a Icexp99]; Id_a = [Id_a ];   Vcm_a = [Vcm_a Vmexp99]; 
Vce_a = [Vce_a Veexp99]; Vc0_a = [Vc0_a Vc0exp99]; Vd0_a = [Vd0_a Vd0exp99]; Tcm_a = [Tcm_a Tmexp99];

vidx = 1:length(Qc_a);
figure(fg0)
plot(Mcyc(j,nodd(vidx)), Qc_a/1e3,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;
plot(Mcyc(j,neven(vidx)),Qd_a/1e3,'og','MarkerFaceColor','g','markersize',6);hold on;axis tight;%title('Capacity fade')
figure(fg0+2)
plot(Mcyc(j,nodd(2:vidx(end))),tc_a(2:end)/3600,'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
figure(fg0+3)
plot(Mcyc(j,nodd(2:vidx(end))),Icm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Average I')
plot(Mcyc(j,nodd(2:vidx(end))),Ic_a(2:end),'om','markersize',6);hold on;axis tight;%title('Average I')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+4)
plot(Mcyc(j,nodd(2:vidx(end))),Vcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(2:vidx(end))),Vce_a(2:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+5)
plot(Mcyc(j,nodd(2:vidx(end))),Tcm_a(2:end),'or','markersize',6);hold on;axis tight;%title('Max T')
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[20 70],'--','color',[1 0.5 0],'linewidth',2);hold on;
figure(fg0+6)
plot(Mcyc(j,nodd(1:vidx(end))),Vc0_a(1:end),'or','markersize',6);hold on;axis tight;%title('Max/End V')
plot(Mcyc(j,nodd(1:vidx(end))),Vd0_a(1:end),'om','markersize',6);hold on;axis tight
% plot((Mcyc(j,nodd(vidx(end)))+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
% break
figure(fg0+7)
plot(Mcyc(j,nodd(2:vidx(end))),(Qc_a(2:end)/1e3)./(tc_a(2:end)/3600),'or','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
vV0 = [vV0 [Vc0_a(1:end); Vd0_a(1:end)]];
% break
count = str2double(expP0(1,:));
fname_p = sprintf('Test%d.csv',count);
[Qc_p,Qd_p, Vdc_p,Vdd_p, tc_p,td_p] = importingPEC_PT(fname_p);
Pc_p = Vdc_p * 5.87/2;
Pd_p = Vdd_p * 5.87/2;

figure(fg0)
plot(Mcyc(j,end-1),Qc_p/1e3,   'or','linewidth',2,'markersize',12);hold on;axis tight;
plot(Mcyc(j,end),  Qd_p(2)/1e3,'og','linewidth',2,'markersize',12);hold on;axis tight;
% plot((Mcyc(j,end)+0.5)*ones(1,2),[0 5],'--','color',[1 0.5 0],'linewidth',2);hold on;
%%
break







for j = 1:n_exp
    for j2 = 1:size(expM0,1);
        if spec0(j2) == 0
        else
%             Qc_a = ;  Qd_a = ;  Vdc_a = ;  Vdd_a = ;  tc_a = ;  
%             td_a = ;  Ic_a = ;  Id_a = ;   Vcm_a = ;  Vce_a = ;  Tcm_a = ;
        end
    %     break
    %     Qc_a_tot = [];
    %     Qd_a_tot = [];
    end


    
    
    
    if j<=size(expP,1)

    end
    
    
    
    
%     count = count+2;
%     
%     if j==n_exp
%         %% Last cycle before Long Performance Test
%         fname_b = 'Test1136.csv';
% %         [Qc_b,Qd_b, Vdc_b,Vdd_b] = importingPEC_ageing(fname_b);
%         [Qc_b,Qd_b, Vdc_b,Vdd_b, tc_b,td_b] = importingPEC_lastC(fname_b, 2);
%         Pc_b = Vdc_b * 5.87;
%         Pd_b = Vdd_b * 5.87/2;
% %         break
%         
% %         %original Cap
%         figure(80)
%         plot(Mcyc(j,end)+[1 3],Qc_b,'or',Mcyc(j,end)+[2 4],Qd_b,'og','markersize',6);hold on;axis tight;title('Capacity fade')
% %         %normalized Cap
% %         figure(80)
% %         plot(Mcyc(j,end)+[1 3],Qc_b/2900/0.9327,'or',Mcyc(j,end)+[2 4],Qd_b/2900/0.9327,'og','markersize',6);hold on;axis tight;title('Capacity fade')
%         figure(81)
%         plot(Mcyc(j,end)+[1 3],Pc_b,'or',Mcyc(j,end)+[2 4],Pd_b,'og','markersize',6);hold on;axis tight;title('Power fade')
% %         figure(82)
% %         plot(Mcyc(j,end)+[1 3],tc_b,'or',Mcyc(j,end)+[2 4],td_b,'og','markersize',6);hold on;axis tight;
% 
%         %% Long Performance Test
% %         fname_a = sprintf('Test11%d.csv',count);
% %         fname_l = 'Test1180.csv';
%         fname_l = 'Test1276.csv';
%         [Qc_l,Qd_l, Vdc_l,Vdd_l, tc_l,td_l] = importingPEC_longPbb(fname_l);
%         Pc_l = Vdc_l * 5.87;
%         Pd_l = Vdd_l * 5.87/2;
% 
% %         %original Cap
%         figure(80)
%         plot(Mcyc(j,end)+4+[1 3],Qc_l,'pr',Mcyc(j,end)+4+[2 4],Qd_l,'pg','markersize',24);hold on;axis tight;title('Capacity fade')
%         %normalized Cap
% %         figure(80)
% %         plot(Mcyc(j,end)+4+[1 3],Qc_l/2900/0.9327,'pr',Mcyc(j,end)+4+[2 4],Qd_l/2900/0.9327,'pg','markersize',24);hold on;axis tight;title('Capacity fade')
%         figure(81)
% %     plot(Mcyc(j,nodd),Vdc_a,'or',Mcyc(j,neven),Vdd_a,'og','markersize',24);hold on;axis tight;title('Power fade')
%         plot(Mcyc(j,end)+4+[1 3],Pc_l,'pr',Mcyc(j,end)+4+[2 4],Pd_l,'pg','markersize',24);hold on;axis tight;title('Power fade')
% %         figure(82)
% %         plot(Mcyc(j,end)+4+[1 3],tc_l,'pr',Mcyc(j,end)+4+[2 4],td_l,'pg','markersize',24);hold on;axis tight;
%     end
%     %in numbers: Charge 1169.953 	1169.943
%                 %Discharge 1169.949	1169.947
end
break

%%
% tag_name = 1;
tag_name = 0;
figure(fg0)
    set(gca,'fontsize',16)
    ylabel('Capacity [Ah]','fontsize',16);xlabel('Half cycle number','fontsize',16)
if tag_name == 1
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
x1=[];  y1=[];
for j = 1:size(xdata,1)
    x0 = xdata{j};  y0 = ydata{j};
    x1 = [x1 x0];   y1 = [y1 y0];
end
% vva = [Mcyc(j,nodd) Mcyc(j,neven)];
% L = strsplit(sprintf('%s\n',repmat('rg,',1,length(x1))),',');  % Letter Labels
L = strsplit(sprintf('%s\n',repmat('s,',1,length(x1))),',');  % Letter Labels
text(x1,y1, L(1:length(x1)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
end
xbar1 = 24.5:24:Mcyc(end,end);  xbar2 = [xbar1' xbar1']';
ybar = repmat([0 5],1,length(xbar2));
% plot(xbar2(:)',ybar,'--','color',[1 0.5 0],'linewidth',2);hold on;
plot((reshape(xbar2(:)',2,length(ybar)/2)')',(reshape(ybar,2,length(ybar)/2)')','--','color',[1 0.5 0],'linewidth',2);hold on;
% break
figure(fg0+2)
    set(gca,'fontsize',16)
    ylabel('Charge Time [h]','fontsize',16);xlabel('Half cycle number','fontsize',16)
if tag_name == 1
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
x1=[];  y1=[];
for j = 1:size(xdata,1)
    x0 = xdata{j};  y0 = ydata{j};
    x1 = [x1 x0];   y1 = [y1 y0];
end
% vva = [Mcyc(j,nodd) Mcyc(j,neven)];
% L = strsplit(sprintf('%s\n',repmat('rg,',1,length(x1))),',');  % Letter Labels
L = strsplit(sprintf('%s\n',repmat('s,',1,length(x1))),',');  % Letter Labels
text(x1,y1, L(1:length(x1)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
end
xbar1 = 24.5:24:Mcyc(end,end);  xbar2 = [xbar1' xbar1']';
ybar = repmat([0 5],1,length(xbar2));
plot((reshape(xbar2(:)',2,length(ybar)/2)')',(reshape(ybar,2,length(ybar)/2)')','--','color',[1 0.5 0],'linewidth',2);hold on;
% break
figure(fg0+3)
    set(gca,'fontsize',16)
    ylabel('Average C-rate [h^{-1}]','fontsize',16);xlabel('Half cycle number','fontsize',16)
if tag_name == 1
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
x1=[];  y1=[];
for j = 1:size(xdata,1)
    x0 = xdata{j};  y0 = ydata{j};
    x1 = [x1 x0];   y1 = [y1 y0];
end
% vva = [Mcyc(j,nodd) Mcyc(j,neven)];
% L = strsplit(sprintf('%s\n',repmat('rg,',1,length(x1))),',');  % Letter Labels
L = strsplit(sprintf('%s\n',repmat('s,',1,length(x1))),',');  % Letter Labels
text(x1,y1, L(1:length(x1)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
end
xbar1 = 24.5:24:Mcyc(end,end);  xbar2 = [xbar1' xbar1']';
ybar = repmat([0 5],1,length(xbar2));
plot((reshape(xbar2(:)',2,length(ybar)/2)')',(reshape(ybar,2,length(ybar)/2)')','--','color',[1 0.5 0],'linewidth',2);hold on;
% break
figure(fg0+4)
    set(gca,'fontsize',16)
    ylabel('Voltage [V]','fontsize',16);xlabel('Half cycle number','fontsize',16)
if tag_name == 1
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
x1=[];  y1=[];
for j = 1:size(xdata,1)
    x0 = xdata{j};  y0 = ydata{j};
    x1 = [x1 x0];   y1 = [y1 y0];
end
% vva = [Mcyc(j,nodd) Mcyc(j,neven)];
% L = strsplit(sprintf('%s\n',repmat('rg,',1,length(x1))),',');  % Letter Labels
L = strsplit(sprintf('%s\n',repmat('s,',1,length(x1))),',');  % Letter Labels
text(x1,y1, L(1:length(x1)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
end
xbar1 = 24.5:24:Mcyc(end,end);  xbar2 = [xbar1' xbar1']';
ybar = repmat([0 5],1,length(xbar2));
plot((reshape(xbar2(:)',2,length(ybar)/2)')',(reshape(ybar,2,length(ybar)/2)')','--','color',[1 0.5 0],'linewidth',2);hold on;
% break
figure(fg0+5)
    set(gca,'fontsize',16)
    ylabel('Max Temperature [dC]','fontsize',16);xlabel('Half cycle number','fontsize',16)
if tag_name == 1
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
x1=[];  y1=[];
for j = 1:size(xdata,1)
    x0 = xdata{j};  y0 = ydata{j};
    x1 = [x1 x0];   y1 = [y1 y0];
end
% vva = [Mcyc(j,nodd) Mcyc(j,neven)];
% L = strsplit(sprintf('%s\n',repmat('rg,',1,length(x1))),',');  % Letter Labels
L = strsplit(sprintf('%s\n',repmat('s,',1,length(x1))),',');  % Letter Labels
text(x1,y1, L(1:length(x1)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
end
xbar1 = 24.5:24:Mcyc(end,end);  xbar2 = [xbar1' xbar1']';
ybar = repmat([20 70],1,length(xbar2));
plot((reshape(xbar2(:)',2,length(ybar)/2)')',(reshape(ybar,2,length(ybar)/2)')','--','color',[1 0.5 0],'linewidth',2);hold on;
% break
figure(fg0+6)
    set(gca,'fontsize',16)
    ylabel('Voltage [V]','fontsize',16);xlabel('Half cycle number','fontsize',16)
if tag_name == 1
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
x1=[];  y1=[];
for j = 1:size(xdata,1)
    x0 = xdata{j};  y0 = ydata{j};
    x1 = [x1 x0];   y1 = [y1 y0];
end
% vva = [Mcyc(j,nodd) Mcyc(j,neven)];
% L = strsplit(sprintf('%s\n',repmat('rg,',1,length(x1))),',');  % Letter Labels
L = strsplit(sprintf('%s\n',repmat('s,',1,length(x1))),',');  % Letter Labels
text(x1,y1, L(1:length(x1)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
end
xbar1 = 24.5:24:Mcyc(end,end);  xbar2 = [xbar1' xbar1']';
ybar = repmat([0 5],1,length(xbar2));
plot((reshape(xbar2(:)',2,length(ybar)/2)')',(reshape(ybar,2,length(ybar)/2)')','--','color',[1 0.5 0],'linewidth',2);hold on;
%
figure(fg0+7)
    set(gca,'fontsize',16)
    ylabel('Capacity per charge time [Ah.h^{-1}]','fontsize',16);xlabel('Half cycle number','fontsize',16)
if tag_name == 1
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
x1=[];  y1=[];
for j = 1:size(xdata,1)
    x0 = xdata{j};  y0 = ydata{j};
    x1 = [x1 x0];   y1 = [y1 y0];
end
% vva = [Mcyc(j,nodd) Mcyc(j,neven)];
% L = strsplit(sprintf('%s\n',repmat('rg,',1,length(x1))),',');  % Letter Labels
L = strsplit(sprintf('%s\n',repmat('s,',1,length(x1))),',');  % Letter Labels
text(x1,y1, L(1:length(x1)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
end
xbar1 = 24.5:24:Mcyc(end,end);  xbar2 = [xbar1' xbar1']';
ybar = repmat([0 5],1,length(xbar2));
% plot(xbar2(:)',ybar,'--','color',[1 0.5 0],'linewidth',2);hold on;
plot((reshape(xbar2(:)',2,length(ybar)/2)')',(reshape(ybar,2,length(ybar)/2)')','--','color',[1 0.5 0],'linewidth',2);hold on;
% break


% cn = 1:N;
% figure(80)
% % subplot(211)
% plot(cn,Qc,'or',cn,Qd,'og','markersize',24);hold on;axis tight;
% % subplot(212)
% % plot(t,V,'linewidth',2);hold on;axis tight;
