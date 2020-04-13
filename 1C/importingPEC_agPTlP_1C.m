clear all
clc

exp1 = '1396';
exp2 = '1397';
exp3 = '1422';
exp4 = '1423';
exp5 = '1424';
exp6 = '1439';
exp7 = '1440';
exp8 = '1441';
exp9 = '1442';
exp10= '1443';

exp11= '1472';
exp12= '1473';
exp13= '1474';
exp14= '1499';
exp15= '1500';
exp16= '1501';
exp17= '1502';
exp18= '1503';
exp19= '1504';
exp20= '1505';

exp21= '1548';
exp22= '1561';

expM = [exp1; exp3; exp5; exp7; exp9;  exp11; exp13; exp15; exp17; exp19; exp21];
expP = [exp2; exp4; exp6; exp8; exp10; exp12; exp14; exp16; exp18; exp20; exp22];
expN = [11; 11; 11; 11; 11; 11; 11; 11; 11; 11; 11];
n_exp = 11;

ncyc_exp = 11 + 1;
vcyc = 1:2*ncyc_exp*n_exp;
Mcyc = reshape(vcyc,2*ncyc_exp,n_exp)';

neven = 2:2:2*ncyc_exp-1;
nodd = 1:2:2*ncyc_exp-2;

fg0 = 60;
QQ = 0;     Q1C_3 = [];     R1C_3 = [];
vvii = [];  idx = 1;
xdata_1c = [];  ydata_1c = [];
Mcycv = [];     Tcm_av = [];    ts = 1;
Mcycv1 = [];    Vc0_av1 = [];
Mcycv2 = [];    Vd0_av2 = [];

dI = 2e-3;    dV = 2e-3;    
MdCpc=0;    MdCpd=0;    MdPc=0;    MdPd=0;    Mdtc=0;
c0 = 0;    g = 1;   

CycCapv = [];   CycTimv = [];   RPTCapv = [];   RPTTimv = [];
CycResv = [];   RPTResv = [];   CycCapChTv = [];
CycChTv = [];   CycChTTimv = [];
CycCapv_b = [];     CycChTv_b = [];
CycCapv_Q1 = [];    RPTCapv_Q1 = [];
CycCapv_Q2 = [];    RPTCapv_Q2 = [];
for j = 1:n_exp
%%
    count = str2double(expM(j,:));
    fname_a = sprintf('Test%d.csv',count);
    [Qc_a,Qd_a, Vdc_a,Vdd_a,Vdd2_a, tc_a,td_a, Icm_a,Ic_a,Id_a,Vcm_a,Vce_a,Vc0_a,Vd0_a,Tcm_a, Qci_a,tci_a] = importingPEC_ageing(fname_a, expN(j));
    Pc_a = Vdc_a * 5.87;    Pd_a = Vdd_a * 5.87/2;
    Pc_a2 = Vdc_a/2.75;     Pd_a2 = Vdd_a/2.75;     Pd2_a2 = Vdd2_a/2.75;

%%
%     %original Cap
    figure(fg0)
    plot(Mcyc(j,nodd(2:ts:end)), (c0+g*Qc_a(2:ts:end))/1e3,'^r','MarkerFaceColor','r','markersize',6);hold on;axis tight;
        plot(Mcyc(j,nodd(2:ts:end)), (c0+g*Qci_a(2:ts:end))/1e3,'^r','MarkerFaceColor','r','MarkerEdgeColor','w','markersize',6);hold on;axis tight;
      dCpc=dI*tc_a(2:ts:end);  %plot(Mcyc(j,nodd(2:ts:end)), (c0+g*Qc_a(2:ts:end)+dCpc)/1e3,'+b','markersize',4);hold on;axis tight;
                               %plot(Mcyc(j,nodd(2:ts:end)), (c0+g*Qc_a(2:ts:end)-dCpc)/1e3,'+b','markersize',4);hold on;axis tight;
      MdCpc=max([MdCpc,dCpc]);
    plot(Mcyc(j,neven(2:ts:end)),(c0+g*Qd_a(2:ts:end))/1e3,'^g','MarkerFaceColor','g','markersize',6);hold on;axis tight;%title('Capacity fade')
      dCpd=dI*td_a(2:ts:end);  %plot(Mcyc(j,neven(2:ts:end)), (c0+g*Qd_a(2:ts:end)+dCpd)/1e3,'+b','markersize',4);hold on;axis tight;
                               %plot(Mcyc(j,neven(2:ts:end)), (c0+g*Qd_a(2:ts:end)-dCpd)/1e3,'+b','markersize',4);hold on;axis tight;
      MdCpd=max([MdCpd,dCpd]);
      CycCap1 = [Qc_a(2:ts:end); Qd_a(2:ts:end)]/1e3;  CycCap2 = CycCap1(:)';  CycCapv = [CycCapv CycCap2];  CycTimv = [CycTimv Mcyc(j,1+2:end-2)];
      CycCapv_b = [CycCapv_b Qci_a(2:ts:end)/1e3];
      
    figure(fg0+1)
    plot(Mcyc(j,nodd(2:ts:end)), Pc_a2(2:ts:end),'^r','MarkerFaceColor','r','markersize',6);hold on;axis tight;
      dPc=(2*dV./Vdc_a(2:ts:end)+2*dI/(2*2.75)).*Pc_a2(2:ts:end);  %plot(Mcyc(j,nodd(2:ts:end)), Pc_a2(2:ts:end)+dPc,'+b','markersize',4);hold on;axis tight;
                                                                   %plot(Mcyc(j,nodd(2:ts:end)), Pc_a2(2:ts:end)-dPc,'+b','markersize',4);hold on;axis tight;
    MdPc=max([MdPc,dPc]);
    plot(Mcyc(j,neven(2:ts:end)),Pd_a2(2:ts:end),'^g','MarkerFaceColor','g','markersize',6);hold on;axis tight;%title('Power fade')
      dPd=(2*dV./Vdd_a(2:ts:end)+2*dI/2.75).*Pd_a2(2:ts:end);  %plot(Mcyc(j,neven(2:ts:end)), Pd_a2(2:ts:end)+dPd,'+b','markersize',4);hold on;axis tight;
                                                               %plot(Mcyc(j,neven(2:ts:end)), Pd_a2(2:ts:end)-dPd,'+b','markersize',4);hold on;axis tight;
    MdPd=max([MdPd,dPd]);
    CycRes1 = [Pc_a2(2:ts:end); Pd_a2(2:ts:end)];  CycRes2 = CycRes1(:)';  CycResv = [CycResv CycRes2];
xdata_1c = [xdata_1c Mcyc(j,nodd(2:end))];
ydata_1c = [ydata_1c tc_a(2:end)];
    
    figure(fg0+2)
    plot(Mcyc(j,nodd(2:ts:end)),tc_a(2:ts:end)/3600,'^r','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
        plot(Mcyc(j,nodd(2:ts:end)),tci_a(2:ts:end)/3600,'^r','MarkerFaceColor','r','MarkerEdgeColor','w','markersize',6);hold on;axis tight;%title('Charge time')
    dtc=1;  %plot(Mcyc(j,nodd(2:ts:end)), (tc_a(2:ts:end)+dtc)/3600,'+b','markersize',4);hold on;axis tight;
            %plot(Mcyc(j,nodd(2:ts:end)), (tc_a(2:ts:end)-dtc)/3600,'+b','markersize',4);hold on;axis tight;
    Mdtc=max([Mdtc,dtc]);
    CycChT1 = tc_a(2:ts:end)/3600;  CycChT2 = CycChT1(:)';  CycChTv = [CycChTv CycChT2];  CycChTTimv = [CycChTTimv Mcyc(j,nodd(1+1:ts:end))];
    CycChTv_b = [CycChTv_b tci_a(2:ts:end)/3600];
    xlabel('Cycle number');ylabel('Time [h]')
    
%     figure(fg0+3)
%     plot(Mcyc(j,nodd(2:ts:end)),Icm_a(2:ts:end),'^r','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Average I')
%     plot(Mcyc(j,nodd(2:ts:end)),Ic_a(2:ts:end),'^r','linewidth',2,'markersize',6);hold on;axis tight;%title('Average I')
    
%     figure(fg0+4)
%     plot(Mcyc(j,nodd(2:ts:end)),Vcm_a(2:ts:end),'^r','linewidth',2,'markersize',6);hold on;axis tight;%title('Max/End V')
%     plot(Mcyc(j,nodd(2:ts:end)),Vce_a(2:ts:end),'^r','MarkerFaceColor','r','markersize',6);hold on;axis tight
    
%     figure(fg0+5)
%     plot(Mcyc(j,nodd(2:ts:end)),Tcm_a(2:ts:end),'^r','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Max T')
%     Mcycv = [Mcycv Mcyc(j,nodd(2:ts:end))];    Tcm_av = [Tcm_av Tcm_a(2:ts:end)];
    
%     figure(fg0+6)
%     plot(Mcyc(j,nodd(1:ts:end)),Vc0_a(1:ts:end),'^r','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Max/End V')
%     plot(Mcyc(j,nodd(1:ts:end)),Vd0_a(1:ts:end),'^r','linewidth',2,'markersize',6);hold on;axis tight
    Mcycv1 = [Mcycv1 Mcyc(j,nodd(1:ts:end))];      Vc0_av1 = [Vc0_av1 Vc0_a(1:ts:end)];
    Mcycv2 = [Mcycv2 Mcyc(j,nodd(1:ts:end))];      Vd0_av2 = [Vd0_av2 Vd0_a(1:ts:end)];
    
%     figure(fg0+7)
%     plot(Mcyc(j,nodd(2:end)),(Qc_a(2:end)/1e3)./(tc_a(2:end)/3600),'^r','MarkerFaceColor','r','markersize',6);hold on;axis tight;%title('Charge time')
%     xlabel('Cycle number');ylabel('Capacity per charge time [Ah.h^{-1}]')
    Qx = [Qc_a'/1e3 Qd_a'/1e3]; Qy = Qx';   Qz=Qy(:);   QQ=QQ(end)+cumsum(Qz);
    CycCapChT1 = (Qc_a(2:end)/1e3)./(tc_a(2:end)/3600);  CycCapChT2 = CycCapChT1(:)';  CycCapChTv = [CycCapChTv CycCapChT2];
    
%     figure(fg0+8)
%     plot(QQ(1:2:length(QQ)),Qc_a/1e3,'^r','MarkerFaceColor','r','markersize',6);hold on;axis tight;
%     plot(QQ(2:2:length(QQ)),Qd_a/1e3,'^g','MarkerFaceColor','g','markersize',6);hold on;axis tight;%title('Capacity fade')
        CycCap1_Q1 = [Qc_a; Qd_a]/1e3;  CycCap2_Q1 = CycCap1_Q1(:)';  CycCapv_Q1 = [CycCapv_Q1 CycCap2_Q1];
        CycCapv_Q2 = [CycCapv_Q2 QQ'];
    
    count = str2double(expP(j,:));
    fname_p = sprintf('Test%d.csv',count);
    [Qc_p,Qd_p, Idc_p,Vdc_p,Idd_p,Vdd_p, tc_p,td_p] = importingPEC_PT(fname_p);

%     %original Cap
    figure(fg0)
    plot(Mcyc(j,end-1),(c0+g*Qc_p)   /1e3,'^r','linewidth',2,'markersize',12);hold on;axis tight;
      dCpc_p=dI*tc_p;  %plot(Mcyc(j,end-1), (c0+g*Qc_p+dCpc_p)/1e3,'+b','markersize',4);hold on;axis tight;
                       %plot(Mcyc(j,end-1), (c0+g*Qc_p-dCpc_p)/1e3,'+b','markersize',4);hold on;axis tight;
      MdCpc=max([MdCpc,dCpc_p]);
    plot(Mcyc(j,end),  (c0+g*Qd_p(2))/1e3,'^g','linewidth',2,'markersize',12);hold on;axis tight;
      dCpd_p=dI*td_p(2);  %plot(Mcyc(j,end), (c0+g*Qd_p(2)+dCpd_p)/1e3,'+b','markersize',4);hold on;axis tight;
                          %plot(Mcyc(j,end), (c0+g*Qd_p(2)-dCpd_p)/1e3,'+b','markersize',4);hold on;axis tight;
      MdCpd=max([MdCpd,dCpd_p]);
    Q1C_1 = [Qc_p; Qd_p(2)];   Q1C_2 = Q1C_1(:)';   Q1C_3 = [Q1C_3 Q1C_2];
    RPTCap1 = [Qc_p; Qd_p(2)]/1e3;  RPTCap2 = RPTCap1(:)';  RPTCapv = [RPTCapv RPTCap2];  RPTTimv = [RPTTimv Mcyc(j,end-1:end)];

    figure(fg0+1)
    plot(Mcyc(j,end-1),Vdc_p/Idc_p,'^r','linewidth',2,'markersize',12);hold on;axis tight;
    dPc_p=(2*dV./Vdc_p+2*dI/Idc_p).*Vdc_p/Idc_p;  %plot(Mcyc(j,end-1), Vdc_p/Idc_p+dPc_p,'+b','markersize',4);hold on;axis tight;
                                                  %plot(Mcyc(j,end-1), Vdc_p/Idc_p-dPc_p,'+b','markersize',4);hold on;axis tight;
    MdPc=max([MdPc,dPc_p]);
    plot(Mcyc(j,end),  Vdd_p(2)/Idd_p(2),'^g','linewidth',2,'markersize',12);hold on;axis tight;%title('Power fade')
    dPd_p=(2*dV./Vdd_p(2)+2*dI/Idd_p(2)).*Vdd_p(2)/Idd_p(2);  %plot(Mcyc(j,end), Vdd_p(2)/Idd_p(2)+dPd_p,'+b','markersize',4);hold on;axis tight;
                                                              %plot(Mcyc(j,end), Vdd_p(2)/Idd_p(2)-dPd_p,'+b','markersize',4);hold on;axis tight;
    MdPd=max([MdPd,dPd_p]);
    R1C_1 = [Vdc_p/Idc_p; Vdd_p(2)/Idd_p(2)];   R1C_2 = R1C_1(:)';   R1C_3 = [R1C_3 R1C_2];
    Qx = [Qc_p'/1e3 Qd_p(2)'/1e3]; Qy = Qx';   Qz=Qy(:);   QQ=QQ(end)+cumsum(Qz);
    RPTRes1 = [Vdc_p/Idc_p; Vdd_p(2)/Idd_p(2)];  RPTRes2 = RPTRes1(:)';  RPTResv = [RPTResv RPTRes2];
    
%     figure(fg0+8)
%     plot(QQ(1:2:length(QQ)),Qc_p/1e3,'^r','linewidth',2,'markersize',12);hold on;axis tight;
%     plot(QQ(2:2:length(QQ)),Qd_p(2)/1e3,'^g','linewidth',2,'markersize',12);hold on;axis tight;%title('Capacity fade')
        RPTCap1_Q1 = [Qc_p; Qd_p(2)]/1e3;  RPTCap2_Q1 = RPTCap1_Q1(:)';  RPTCapv_Q1 = [RPTCapv_Q1 RPTCap2_Q1];
        RPTCapv_Q2 = [RPTCapv_Q2 QQ'];
        
%     figure(fg0+9)
%     plot(idx,Qc_p   /1e3,'^r','linewidth',2,'markersize',12);hold on;axis tight;  idx = idx+1;
%     plot(idx,Qd_p(2)/1e3,'^g','linewidth',2,'markersize',12);hold on;axis tight;  idx = idx+1;
    vvii = [vvii; Mcyc(j,end-1); Mcyc(j,end)];
% pause

    if j == 5
        %% Long Performance Test
        fname_l = 'Test1444.csv';
        [Qc_l,Qd_l, Vdc_l,Vdd_l, tc_l,td_l] = importingPEC_longPbb(fname_l);
        Pc_l = Vdc_l * 5.87;
        Pd_l = Vdd_l * 5.87/2;

%         %original Cap
        figure(fg0)
        plot(Mcyc(j,end)+[1 3],(c0+g*Qc_l)/1e3,'pr',Mcyc(j,end)+[2 4],(c0+g*Qd_l)/1e3,'pg','markersize',24);hold on;axis tight;%title('Capacity fade')
    elseif j == 10
        %% Long Performance Test
        fname_l = 'Test1531.csv';
        [Qc_l,Qd_l, Vdc_l,Vdd_l, tc_l,td_l] = importingPEC_longPbb(fname_l);
        Pc_l = Vdc_l * 5.87;
        Pd_l = Vdd_l * 5.87/2;

%         %original Cap
        figure(fg0)
        plot(Mcyc(j,end)+[1 3],(c0+g*Qc_l)/1e3,'pr',Mcyc(j,end)+[2 4],(c0+g*Qd_l)/1e3,'pg','markersize',24);hold on;axis tight;%title('Capacity fade')
    end

end
num2str([MdCpc    MdCpd    MdPc    MdPd    Mdtc],'%2.5e\t')

save QR_1C_e.mat Q1C_3 R1C_3

results1C.xdata.cyc = CycTimv;
results1C.xdata.rpt = RPTTimv;
results1C.xdata.chg = CycChTTimv;
results1C.ydata.capacity_cyc = CycCapv;
results1C.ydata.capacity_cycb = CycCapv_b;
results1C.ydata.capacity_rpt = RPTCapv;
results1C.ydata.resistance_cyc = CycResv;
results1C.ydata.resistance_rpt = RPTResv;
results1C.ydata.chgtime = CycChTv;
results1C.ydata.chgtimeb = CycChTv_b;
results1C.ydata.capacity_chgtime = CycCapChTv;

results1C.ydata.capacity_cyc_thput = CycCapv_Q1;
results1C.ydata.capacity_rpt_thput = RPTCapv_Q1;
results1C.xdata.cyc_thput = CycCapv_Q2;
results1C.xdata.rpt_thput = RPTCapv_Q2;

save data1C_e.mat results1C

%%
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
L = strsplit(sprintf('%s\n',repmat('s,',1,length(x1))),',');  % Letter Labels
text(x1,y1, L(1:length(x1)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
end
xbar1 = 24.5:24:Mcyc(end,end);  xbar2 = [xbar1' xbar1']';
ybar = repmat([0 5],1,length(xbar2));
plot((reshape(xbar2(:)',2,length(ybar)/2)')',(reshape(ybar,2,length(ybar)/2)')','--','color',[1 0.5 0],'linewidth',2);hold on;

figure(fg0+1)
    set(gca,'fontsize',16)
    ylabel('Resistance []','fontsize',16);xlabel('Half cycle number','fontsize',16)
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
L = strsplit(sprintf('%s\n',repmat('s,',1,length(x1))),',');  % Letter Labels
text(x1,y1, L(1:length(x1)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
end
xbar1 = 24.5:24:Mcyc(end,end);  xbar2 = [xbar1' xbar1']';
ybar = repmat([0 5],1,length(xbar2));
plot((reshape(xbar2(:)',2,length(ybar)/2)')',(reshape(ybar,2,length(ybar)/2)')','--','color',[1 0.5 0],'linewidth',2);hold on;

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
L = strsplit(sprintf('%s\n',repmat('s,',1,length(x1))),',');  % Letter Labels
text(x1,y1, L(1:length(x1)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
end
xbar1 = 24.5:24:Mcyc(end,end);  xbar2 = [xbar1' xbar1']';
ybar = repmat([0 5],1,length(xbar2));
plot((reshape(xbar2(:)',2,length(ybar)/2)')',(reshape(ybar,2,length(ybar)/2)')','--','color',[1 0.5 0],'linewidth',2);hold on;
% break

% figure(fg0+3)
%     set(gca,'fontsize',16)
%     ylabel('Average C-rate [h^{-1}]','fontsize',16);xlabel('Half cycle number','fontsize',16)
% if tag_name == 1
% h = gcf; %current figure handle
% axesObjs = get(h, 'Children');  %axes handles
% dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
% xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
% ydata = get(dataObjs, 'YData');
% x1=[];  y1=[];
% for j = 1:size(xdata,1)
%     x0 = xdata{j};  y0 = ydata{j};
%     x1 = [x1 x0];   y1 = [y1 y0];
% end
% L = strsplit(sprintf('%s\n',repmat('s,',1,length(x1))),',');  % Letter Labels
% text(x1,y1, L(1:length(x1)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
% end
% xbar1 = 24.5:24:Mcyc(end,end);  xbar2 = [xbar1' xbar1']';
% ybar = repmat([-5 0],1,length(xbar2));
% plot((reshape(xbar2(:)',2,length(ybar)/2)')',(reshape(ybar,2,length(ybar)/2)')','--','color',[1 0.5 0],'linewidth',2);hold on;

% figure(fg0+4)
%     set(gca,'fontsize',16)
%     ylabel('Voltage [V]','fontsize',16);xlabel('Half cycle number','fontsize',16)
% if tag_name == 1
% h = gcf; %current figure handle
% axesObjs = get(h, 'Children');  %axes handles
% dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
% xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
% ydata = get(dataObjs, 'YData');
% x1=[];  y1=[];
% for j = 1:size(xdata,1)
%     x0 = xdata{j};  y0 = ydata{j};
%     x1 = [x1 x0];   y1 = [y1 y0];
% end
% L = strsplit(sprintf('%s\n',repmat('s,',1,length(x1))),',');  % Letter Labels
% text(x1,y1, L(1:length(x1)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
% end
% xbar1 = 24.5:24:Mcyc(end,end);  xbar2 = [xbar1' xbar1']';
% ybar = repmat([0 5],1,length(xbar2));
% plot((reshape(xbar2(:)',2,length(ybar)/2)')',(reshape(ybar,2,length(ybar)/2)')','--','color',[1 0.5 0],'linewidth',2);hold on;

% figure(fg0+5)
%     set(gca,'fontsize',16)
%     ylabel('Max Temperature [dC]','fontsize',16);xlabel('Half cycle number','fontsize',16)
% if tag_name == 1
% h = gcf; %current figure handle
% axesObjs = get(h, 'Children');  %axes handles
% dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
% xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
% ydata = get(dataObjs, 'YData');
% x1=[];  y1=[];
% for j = 1:size(xdata,1)
%     x0 = xdata{j};  y0 = ydata{j};
%     x1 = [x1 x0];   y1 = [y1 y0];
% end
% L = strsplit(sprintf('%s\n',repmat('s,',1,length(x1))),',');  % Letter Labels
% text(x1,y1, L(1:length(x1)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
% end
% xbar1 = 24.5:24:Mcyc(end,end);  xbar2 = [xbar1' xbar1']';
% ybar = repmat([20 70],1,length(xbar2));
% plot((reshape(xbar2(:)',2,length(ybar)/2)')',(reshape(ybar,2,length(ybar)/2)')','--','color',[1 0.5 0],'linewidth',2);hold on;

% figure(fg0+6)
%     set(gca,'fontsize',16)
%     ylabel('Voltage [V]','fontsize',16);xlabel('Half cycle number','fontsize',16)
% if tag_name == 1
% h = gcf; %current figure handle
% axesObjs = get(h, 'Children');  %axes handles
% dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
% xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
% ydata = get(dataObjs, 'YData');
% x1=[];  y1=[];
% for j = 1:size(xdata,1)
%     x0 = xdata{j};  y0 = ydata{j};
%     x1 = [x1 x0];   y1 = [y1 y0];
% end
% L = strsplit(sprintf('%s\n',repmat('s,',1,length(x1))),',');  % Letter Labels
% text(x1,y1, L(1:length(x1)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
% end
% xbar1 = 24.5:24:Mcyc(end,end);  xbar2 = [xbar1' xbar1']';
% ybar = repmat([0 5],1,length(xbar2));
% plot((reshape(xbar2(:)',2,length(ybar)/2)')',(reshape(ybar,2,length(ybar)/2)')','--','color',[1 0.5 0],'linewidth',2);hold on;

% figure(fg0+7)
%     set(gca,'fontsize',16)
%     ylabel('Capacity per charge time [Ah.h^{-1}]','fontsize',16);xlabel('Half cycle number','fontsize',16)
% if tag_name == 1
% h = gcf; %current figure handle
% axesObjs = get(h, 'Children');  %axes handles
% dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
% xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
% ydata = get(dataObjs, 'YData');
% x1=[];  y1=[];
% for j = 1:size(xdata,1)
%     x0 = xdata{j};  y0 = ydata{j};
%     x1 = [x1 x0];   y1 = [y1 y0];
% end
% L = strsplit(sprintf('%s\n',repmat('s,',1,length(x1))),',');  % Letter Labels
% text(x1,y1, L(1:length(x1)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
% end
% xbar1 = 24.5:24:Mcyc(end,end);  xbar2 = [xbar1' xbar1']';
% ybar = repmat([0 5],1,length(xbar2));
% plot((reshape(xbar2(:)',2,length(ybar)/2)')',(reshape(ybar,2,length(ybar)/2)')','--','color',[1 0.5 0],'linewidth',2);hold on;

% figure(fg0+8)
%     set(gca,'fontsize',16)
%     ylabel('Capacity [Ah]','fontsize',16);xlabel('Throughput [Ah]','fontsize',16)
% if tag_name == 1
% h = gcf; %current figure handle
% axesObjs = get(h, 'Children');  %axes handles
% dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
% xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
% ydata = get(dataObjs, 'YData');
% x1=[];  y1=[];
% for j = 1:size(xdata,1)
%     x0 = xdata{j};  y0 = ydata{j};
%     x1 = [x1 x0];   y1 = [y1 y0];
% end
% L = strsplit(sprintf('%s\n',repmat('s,',1,length(x1))),',');  % Letter Labels
% text(x1,y1, L(1:length(x1)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
% end
% xbar1 = 24.5:24:Mcyc(end,end);  xbar2 = [xbar1' xbar1']';
% ybar = repmat([0 5],1,length(xbar2));
% plot((reshape(xbar2(:)',2,length(ybar)/2)')',(reshape(ybar,2,length(ybar)/2)')','--','color',[1 0.5 0],'linewidth',2);hold on;




