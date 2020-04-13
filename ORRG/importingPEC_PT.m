function [QQc,QQd, VVdc,VVdd, ttc,ttd] = importingPEC_PT(fname)
% fname = 'Test1127.csv';

R1 = 22-1;
% C1 = 8-1;
C1 = 8-1-1;

M = csvread(fname,R1,C1);

st = M(:,1);

t = M(:,2);
V = M(:,3);
% I = M(:,4);

Qc_cum = M(:,5);
Qd_cum = M(:,6);

%discharge 1
std = 1;
[d1,s1,j1] = find(st==std);
wtest = [d1(1):d1(end)]';
d2 = [d1; NaN*ones(length(wtest)-length(d1),1)];
w3 = d2 - wtest;
d3 = d2(w3==0);
Qd = Qd_cum(d3(end));
Vdd = V(d3(1)-1) - V(d3(1));
td = t(d3(end)) - t(d3(1));

%charge
stc = 2;
[c1,r1,i1] = find(st==stc);
vtest = [c1(1):c1(end)]';
c2 = [c1; NaN*ones(length(vtest)-length(c1),1)];
v3 = c2 - vtest;
c3 = c2(v3==0);
Qc = Qc_cum(c3(end));
Vdc = V(c3(1)) - V(c3(1)-1);
tc = t(c3(end)) - t(c3(1));

%discharge 2
stdb = 4;
[d1b,s1b,j1b] = find(st==stdb);
wtestb = [d1b(1):d1b(end)]';
d2b = [d1b; NaN*ones(length(wtestb)-length(d1b),1)];
w3b = d2b - wtestb;
d3b = d2b(w3b==0);
Qdb = Qd_cum(d3b(end))-Qd;
Vddb = V(d3b(1)-1) - V(d3b(1));
tdb = t(d3b(end)) - t(d3b(1));

% N = 3;
% cn = 1:N;
% figure(80)
% % subplot(211)
% plot(cn(1),Qd,'og',cn(2),Qc,'or',cn(3),Qdb,'og','markersize',24);hold on;axis tight;
% % subplot(212)
% % plot(t,V,'linewidth',2);hold on;axis tight;

QQc = Qc;
QQd = [Qd Qdb];

VVdc = Vdc;
VVdd = [Vdd Vddb];

ttc = tc;
ttd = [td tdb];



