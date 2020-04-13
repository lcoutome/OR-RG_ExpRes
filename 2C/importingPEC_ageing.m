function [Qc,Qd, Vdc,Vdd,Vdd2, tc,td, Icm,Ica,Ida,Vcm,Vce,Vc0,Vd0,Tcm, Qci,tci] = importingPEC_ageing(fname, N)
% fname = 'Test1126.csv';

R1 = 22-1;
% C1 = 8-1;
C1 = 8-1-1;

M = csvread(fname,R1,C1);

st = M(:,1);

t = M(:,2);
V = M(:,3);
I = M(:,4);

Qc_cum = M(:,5);
Qd_cum = M(:,6);

T = M(:,7);

stc = 1;
[c1,r1,i1] = find(st==stc);

OneC = -2.75;
% OneC = -2.75;

%charge
% N=11;
Qc = zeros(1,N);
Icm = zeros(1,N);
Vdc = zeros(1,N);
Vc0 = zeros(1,N);
Vcm = zeros(1,N);
Vce = zeros(1,N);
Tcm = zeros(1,N);
tc = zeros(1,N);

Qci = zeros(1,N);
tci = zeros(1,N);

for j=1:N
vtest = [c1(1):c1(end)]';

c2 = [c1; NaN*ones(length(vtest)-length(c1),1)];
v3 = c2 - vtest;
c3 = c2(v3==0);

if j==1
    Qc(j) = Qc_cum(c3(end));
%     tc(j) = t(c3(end));
    tc(j) = t(c3(end)) - t(c3(1));
        Vv = V(c3(1):c3(end));  
        Qq = Qc_cum(c3(1):c3(end));
        tt = t(c3(1):c3(end));
        idn = find(Vv>=(4.2-1e-4));
%         Qci(j) = Qc_cum(idn(1)); tci(j) = t(idn(1)) - t(c3(1));
        Qci(j) = Qq(idn(1)); 
        tci(j) = tt(idn(1)) - tt(1);
elseif j>1
    Qc(j) = Qc_cum(c3(end)) - sum(Qc(1:j-1));
%     tc(j) = t(c3(end)) - sum(tc(1:j-1));
    tc(j) = t(c3(end)) - t(c3(1));
        Vv = V(c3(1):c3(end));  
        Qq = Qc_cum(c3(1):c3(end)) - sum(Qc(1:j-1));
        tt = t(c3(1):c3(end));
        idn = find(Vv>=(4.2-1e-4));
%         Qci(j) = Qc_cum(idn(1)); tci(j) = t(idn(1)) - t(c3(1));
        Qci(j) = Qq(idn(1)); 
        tci(j) = tt(idn(1)) - tt(1);
end
% Vv(idn(1))
% tt(1)
% tt(end)
% tci(j)
% Qci(j)
% pause
Vdc(j) = V(c3(1)) - V(c3(1)-1);

Vc0(j) = V(c3(1)-1);
Icm(j) = max(I(c3(1):c3(end))) / OneC;
Vcm(j) = max(V(c3(1):c3(end)));
Vce(j) = V(c3(end));
Tcm(j) = max(T(c3(1):c3(end)));

%update
% c1 = c1(c3(end)+1:end);
c1 = c1(c3(end)-c1(1)+1+1:end);
% vtest = vtest(c3(end)+1:end);
end
Ica = -(Qc ./ (tc/3600)) / 2900;

%discharge
Qd = zeros(1,N);
Vdd = zeros(1,N);
Vdd2 = zeros(1,N);
Vd0 = zeros(1,N);
Tdm = zeros(1,N);
td = zeros(1,N);
std = 3;
[d1,s1,j1] = find(st==std);
for j=1:N
wtest = [d1(1):d1(end)]';

d2 = [d1; NaN*ones(length(wtest)-length(d1),1)];
w3 = d2 - wtest;
d3 = d2(w3==0);

if j==1
    Qd(j) = Qd_cum(d3(end));
%     td(j) = t(d3(end));
    td(j) = t(d3(end)) - t(d3(1));
elseif j>1
    Qd(j) = Qd_cum(d3(end)) - sum(Qd(1:j-1));
%     td(j) = t(d3(end)) - sum(td(1:j-1));
    td(j) = t(d3(end)) - t(d3(1));
end
Vdd(j) = V(d3(1)-1) - V(d3(1));
Vdd2(j) = V(d3(end)+1) - V(d3(end));

Vd0(j) = V(d3(1)-1);
Tdm(j) = max(T(d3(1):d3(end)));

%update
% c1 = c1(c3(end)+1:end);
d1 = d1(d3(end)-d1(1)+1+1:end);
% vtest = vtest(c3(end)+1:end);
end
Ida = -(Qd ./ (td/3600)) / 2900;

% cn = 1:N;
% figure(80)
% % subplot(211)
% plot(cn,Qc,'or',cn,Qd,'og','markersize',24);hold on;axis tight;
% % subplot(212)
% % plot(t,V,'linewidth',2);hold on;axis tight;

