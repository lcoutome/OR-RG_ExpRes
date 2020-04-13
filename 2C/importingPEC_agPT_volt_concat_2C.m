clear all
clc

% exp1 = '1407';
% % exp2 = '1409';
% exp2 = '1509';

exp1 = '1407';
exp2 = '1409';
exp3 = '1411';
exp4 = '1446';
exp5 = '1457';

exp6 = '1485';
exp7 = '1487';
exp8 = '1509';
exp9 = '1529';

expM1 = exp1;
expM2 = exp2;
expM3 = exp3;
expM4 = exp4;
expM5 = exp5;
expM6 = exp6;
expM7 = exp7;
expM8 = exp8;
expM9 = exp9;
expM_cell = {expM1;expM2;expM3;expM4;expM5;expM6;expM7;expM8;expM9};
co1 = 0;
co2 = 0;
co3 = 0;
co4 = 0;
co5 = 0;
co6 = 0;
co7 = 0;
co8 = 0;
co9 = 0;
co_cell = {co1;co2;co3;co4;co5;co6;co7;co8;co9};

n_exp = size(expM_cell,1);
ncyc_exp = 11 + 1;

%       5      5      5      5      5      5      9       5      6
tti = [48906; 50573; 51440; 70668; 66190; 51521; 102646; 54036; 70570];
ttf = [57655; 59704; 61108; 80664; 76256; 61340; 112612; 63987; 80230];

N=11;
a = rand(N*n_exp+1,3);
i = 1;  l = 1;
idxf = 0;  lh = 2;
tnew = 0;   mI=0;MI=0;  mV=3.7;MV=3.7;  MT=27;
% ba = [ones(9,1) linspace(0,1,9)' zeros(9,1)]; %red-yellow [1 0 0]-[1 1 0]
% bb = [[linspace(0.5,0,5) zeros(1,4)]' ones(9,1) [zeros(1,4) linspace(0,0.5,5)]']; %greens	[0.5 1 0]-[0 1 0]-[0 1 0.5]
bc = [linspace(1,0,9)' zeros(9,1) ones(9,1)]; %pink-blue	[1 0 1]-[0 0 1]

% break
OneC = -2.75;
for l = 1:2;%n_exp;
% l=1;
    expM = expM_cell{l};
    co = co_cell{l};
    for j = 1:length(co);%n_exp
        count = str2double(expM(j,:));
        fname_a = sprintf('Test%d.csv',count);
        [t,V,I,T,Qc_cum,Qd_cum] = importingPEC_ageing_volt_concat(fname_a);
        T = adjs_T(T,t);
%         break
        if co(j)==0
            co2 = length(t);
        else
            co2 = co(j);
        end
        t = t(1:co2); V = V(1:co2); I = I(1:co2); T = T(1:co2);
        Qc_cum = Qc_cum(1:co2);   Qd_cum = Qd_cum(1:co2);
        t = tnew + t;

        figure(70+idxf)
        subplot(311)
        plot([[1:1500] 1500+t']/3600,[ones(1,1500)*I(1) I']/OneC,'--','color',a(i,:),'linewidth',lh);hold on;axis tight;
        mI = min([mI min(I/OneC)]);      MI = max([MI max(I/OneC)]);
        subplot(312)
        plot([[1:1500] 1500+t']/3600,[ones(1,1500)*V(1) V'],'--','color',a(i,:),'linewidth',lh);hold on;axis tight;
        mV = min([mV min(V)]);      MV = max([MV max(V)]);
        subplot(313)
        plot([[1:1500] 1500+t']/3600,[ones(1,1500)*T(1) T'],'--','color',a(i,:),'linewidth',lh);hold on;axis tight;
        MT = max([MT max(T)]);
        
        i = i+1;
        count = count+1;
        tnew = t(end)+1;
    end
    subplot(311);ylim([mI-1/2.9 MI+1/2.9])
        set(gca,'position',[0.13,0.692,0.4,0.28])
        set(gca,'fontsize',16)
%         set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        ylabel('C-rate [h^{-1}]','fontsize',16)
    subplot(312);ylim([mV-0.1 MV+0.1])
        set(gca,'position',[0.13,0.396,0.4,0.28])
        set(gca,'fontsize',16)
        set(gca,'xticklabel',[])
        ylabel('Voltage [V]','fontsize',16)
    subplot(313);ylim([23 MT+1])
        set(gca,'position',[0.13,0.101,0.4,0.28])
        set(gca,'fontsize',16)
        ylabel('Temperature [°C]','fontsize',16);xlabel('Time [h]','fontsize',16)
    tnew = 0;   mI=0;MI=0;  mV=3.7;MV=3.7;  MT=27;

            tex = [[1:1500] 1500+t'];
            Iex = [ones(1,1500)*I(1) I'];
            Vex = [ones(1,1500)*V(1) V'];
            Tex = [ones(1,1500)*T(1) T'];
            tte = tex(tti(l):ttf(l)) - tex(tti(l)) + tex(1);    tte2 = tte(300:end)-tte(300)+tte(1); 
            IIe = Iex(tti(l):ttf(l));                           IIe2 = IIe(300:end); 
            VVe = Vex(tti(l):ttf(l));                           VVe2 = VVe(300:end); 
            TTe = Tex(tti(l):ttf(l));                           TTe2 = TTe(300:end); 
            figure(41)
            subplot(311)
%             plot(tte/3600,IIe/OneC,'--','color',a(i-1,:),'linewidth',lh);hold on;axis tight;   %COLOR CODE FOR CHECKING
            plot(tte2/3600,IIe2/OneC,'--','color',bc(l,:),'linewidth',lh);hold on;axis tight;
            mI = min([mI min(IIe/OneC)]);      MI = max([MI max(IIe/OneC)]);
            subplot(312)
%             plot(tte/3600,VVe,'--','color',a(i-1,:),'linewidth',lh);hold on;axis tight;   %COLOR CODE FOR CHECKING
            plot(tte2/3600,VVe2,'--','color',bc(l,:),'linewidth',lh);hold on;axis tight;
            mV = min([mV min(VVe)]);      MV = max([MV max(VVe)]);
            subplot(313)
%             plot(tte/3600,TTe,'--','color',a(i-1,:),'linewidth',lh);hold on;axis tight;   %COLOR CODE FOR CHECKING
            plot(tte2/3600,TTe2,'--','color',bc(l,:),'linewidth',lh);hold on;axis tight;
            MT = max([MT max(TTe)]);
%             disp(l)
%             disp('done')
%             pause
%             disp('wait')

    idxf = idxf+1;
end

% mI_tot = mean


% cn = 1:N;
% figure(80)
% % subplot(211)
% plot(cn,Qc,'or',cn,Qd,'og','markersize',24);hold on;axis tight;
% % subplot(212)
% % plot(t,V,'linewidth',2);hold on;axis tight;


    subplot(311);ylim([mI-1/2.9 MI+1/2.9])
        set(gca,'position',[0.13,0.692,0.4,0.28])
        set(gca,'fontsize',16)
%         set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        ylabel('C-rate [h^{-1}]','fontsize',16)
    subplot(312);ylim([mV-0.1 MV+0.1])
        set(gca,'position',[0.13,0.396,0.4,0.28])
        set(gca,'fontsize',16)
        set(gca,'xticklabel',[])
        ylabel('Voltage [V]','fontsize',16)
    subplot(313);ylim([23 MT+1])
        set(gca,'position',[0.13,0.101,0.4,0.28])
        set(gca,'fontsize',16)
        ylabel('Temperature [°C]','fontsize',16);xlabel('Time [h]','fontsize',16)
        
        
