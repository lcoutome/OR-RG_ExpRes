function [Qc,Qd, Vdc,Vdd, tc,td] = importingPEC_longPbb(fname)
% fname = 'Test1180.csv';

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

stc = 1;
[c1,r1,i1] = find(st==stc);

%charge
N=2;
Qc = zeros(1,N);
Vdc = zeros(1,N);
tc = zeros(1,N);
for j=1:N
vtest = [c1(1):c1(end)]';

c2 = [c1; NaN*ones(length(vtest)-length(c1),1)];
v3 = c2 - vtest;
c3 = c2(v3==0);

if j==1
    Qc(j) = Qc_cum(c3(end));
    tc(j) = t(c3(end)) - t(c3(1));
elseif j>1
    Qc(j) = Qc_cum(c3(end)) - sum(Qc(1:j-1));
    tc(j) = t(c3(end)) - t(c3(1));
end
Vdc(j) = V(c3(1)) - V(c3(1)-1);

%update
% c1 = c1(c3(end)+1:end);
c1 = c1(c3(end)-c1(1)+1+1:end);
% vtest = vtest(c3(end)+1:end);
end

%discharge
Qd = zeros(1,N);
Vdd = zeros(1,N);
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
    td(j) = t(d3(end)) - t(d3(1));
elseif j>1
    Qd(j) = Qd_cum(d3(end)) - sum(Qd(1:j-1));
    td(j) = t(d3(end)) - t(d3(1));
end
Vdd(j) = V(d3(1)-1) - V(d3(1));

%update
% c1 = c1(c3(end)+1:end);
d1 = d1(d3(end)-d1(1)+1+1:end);
% vtest = vtest(c3(end)+1:end);
end

% cn = 1:N;
% figure(80)
% % subplot(211)
% plot(cn,Qc,'or',cn,Qd,'og','markersize',24);hold on;axis tight;
% % subplot(212)
% % plot(t,V,'linewidth',2);hold on;axis tight;

