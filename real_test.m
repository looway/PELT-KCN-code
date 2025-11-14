clear
close all
load('test3.2_with.mat')

L=length(m);
fs = 10e6;%ADC SampleRate (Hz)
c=3e8;
t=linspace(0,(L-1)/fs*1e6,L);
figure
subplot(2,1,1)
plot(t,m);hold on
% WEN法去噪
%找包络
up_Envelope = WEN(m,5);
plot(t,up_Envelope,'r--');
legend('signal','Envelope');
hold off
xlabel('时间（μs）')
ylabel('幅值')
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))
% 计算变点个数
beta = 0.4;%控制参数
Tc = beta*max(up_Envelope);
tesline=zeros(size(m));
for i=1:length(up_Envelope)
 if up_Envelope(i)>=Tc
    tesline(i)=1;
 else
    tesline(i)=0;
 end
    
end
Y1 = diff(tesline);
Y2 = [0 Y1];
cptn = find(Y2);
subplot(2,1,2)
plot(t,Y2);
xlabel('时间（μs）')
ylabel('幅值')
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))
ipt = findchangepts(m,'Statistic','rms','MaxNumChanges',length(cptn));
% ipt = [325 345];%for3.3
% toc
M = reshape(ipt,2,[]).';
maxl=M(:,1);
minl=M(:,2);
% colum=length(ipt)/2;
% M = reshape(ipt, [2, colum]).';
figure
plot(t,m);
ipt=ipt';
xpos = transpose(double(ipt(1:end))*ones(1,2));
ypos = diag(ylim)*ones(2,numel(ipt));
line(t(xpos),ypos,'Color','green');
xlabel('时间（μs）')
ylabel('幅值')
title('原信号时域图')
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))

Nft=L;
% sweepSlope = 29.9817 * 1e12; %Slope const (MHz/usec)
sweepSlope = 10 * 1e12; %Slope const (MHz/usec)
rangeIdxToMeters = c * fs / (2 * sweepSlope * Nft);
rnggrid = linspace(0, (L/2) * rangeIdxToMeters, L/2+1);%求距离维的坐标标度
% guard = 6;
% train = 30;
% Pf = 1e-2;

%% AWEN
beta = 0.2;%控制参数
Tc = (1+beta)*min(up_Envelope);
wk=zeros(size(m));
for i=1:length(m)
    if abs(m(i))<Tc
        wk(i)=1/Tc;
    else 
        wk(i)=1/up_Envelope(i);
    end
end
AWENsig=m.*wk;
figure
plot(t,AWENsig);
xlabel('时间（μs）')
ylabel('幅值')
title('AWEN归一化时域图');
Y=fft(AWENsig, L);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
% [detected,th] = myCFAR(P1,guard,train,Pf);
% detindx=find(detected);
cfar = phased.CFARDetector('NumTrainingCells',40,'NumGuardCells',2);
exp_pfa = 0.2e-2;
cfar.ThresholdFactor = 'Auto';
cfar.ProbabilityFalseAlarm = exp_pfa;
cfar.ThresholdOutputPort = true;
cfar.Method='SOCA';
P1=P1.';
[detected,th] = cfar(P1,1:length(P1));
targetidx=find(detected);
sigA=sum(P1(targetidx));
noiseA=sum(P1)-sigA;
SNR3=20*log10(sigA/noiseA);
formatSpec = "AWEN距离谱\nSNR is: %.3fdB";
str = sprintf(formatSpec,SNR3);

figure
plot(rnggrid,P1,rnggrid,th,'g--',...
    rnggrid(targetidx),P1(targetidx),'o')
title(str)
xlabel('距离 (米)')
ylabel('幅值')
h=legend('距离谱','CFAR阈值','目标');
set(h, 'Box', 'off')
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))



%% AR
%先测出干扰的起始与结束时间点
% [M,M2]=findjumping(m);
% if length(M)~=length(M2)
%     return
% end
% x=abs(m);
% C = smoothdata(x,'gaussian',100);%使用长度为 40 的较大窗口对原始数据进行平滑处理
% Y1 = diff(C);
% Y2 = [NaN Y1];
% M = find(Y2==max(Y2));
% M2 = find(Y2==min(Y2));

% M = 100;    % Point at which to start predicting
% M2 = 1792;    % Point at which to end predicting
P = length(m);    % Number of samples in the extrapolated time series
% t = 1:P;
m_Original=m;
for k=1:length(minl)
for i=maxl(k):minl(k)
    m(i)=NaN;
end
end
lb = fillgaps(m,3000,200);

% lb=myfilgap(AWENsig,200,M,M2);

figure
subplot(2,1,1)
plot(t,m)
xpos = transpose(double(ipt(1:end))*ones(1,2));
ypos = diag(ylim)*ones(2,numel(ipt));
line(t(xpos),ypos,'Color','green');
xlabel('时间（μs）')
ylabel('幅值')
title('Original');
subplot(2,1,2)
plot(t,lb)
xpos = transpose(double(ipt(1:end))*ones(1,2));
ypos = diag(ylim)*ones(2,numel(ipt));
line(t(xpos),ypos,'Color','green');
xlabel('时间（μs）')
ylabel('幅值')
title('Reconstructed')
% legend('Original','Reconstructed')
Y=fft(lb, L);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
% [detected,th] = myCFAR(P1,guard,train,Pf);
% detindx=find(detected);

cfar = phased.CFARDetector('NumTrainingCells',40,'NumGuardCells',2);
exp_pfa = 1e-3;
cfar.ThresholdFactor = 'Auto';
cfar.ProbabilityFalseAlarm = exp_pfa;
cfar.ThresholdOutputPort = true;
cfar.Method='SOCA';
P1=P1.';
CUTIdx=1:length(P1);
[detected,th] = cfar(P1,1:length(P1));
targetidx=find(detected);
targetidx=targetidx(1);
sigA=sum(P1(targetidx));
noiseA=sum(P1)-sigA;
SNR4=20*log10(sigA/noiseA);
formatSpec = "AR距离谱\nSNR is: %.3fdB";
str = sprintf(formatSpec,SNR4);

figure
plot(rnggrid,P1,rnggrid,th,'g--',...
    rnggrid(targetidx),P1(targetidx),'o')
title(str)
xlabel('距离 (米)')
ylabel('幅值')
h=legend('距离谱','CFAR阈值','目标');
set(h, 'Box', 'off')
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))
%% 
Y=fft(m_Original, L);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
% [detected,th] = myCFAR(P1,guard,train,Pf);
% detindx=find(detected);
sigA=sum(P1(targetidx));
noiseA=sum(P1)-sigA;
SNR1=20*log10(sigA/noiseA);
formatSpec = "无处理距离谱\nSNR is: %.3fdB";
str = sprintf(formatSpec,SNR1);
figure
% plot(rnggrid,P1,rnggrid,th,...
%     rnggrid(targetidx),P1(targetidx),'o')
plot(rnggrid,P1,...
    rnggrid(targetidx),P1(targetidx),'o')
title(str)
xlabel('距离 (米)')
ylabel('幅值')
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))
%% 保存图片
saveas(1,'1.jpg');
saveas(2,'2.jpg');
saveas(3,'3.jpg');
saveas(4,'4.jpg');
saveas(5,'5.jpg');
saveas(6,'6.jpg');
saveas(7,'7.jpg');