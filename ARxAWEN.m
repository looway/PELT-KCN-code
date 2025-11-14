clear
close all
%% 雷达参数设置
% Compute hardware parameters from specified long-range requirements
fc = 77e9;                                  % Center frequency (Hz)
c = physconst('LightSpeed');                % Speed of light in air (m/s)
lambda = c/fc;                              % Wavelength (m)

% Set the chirp duration to be 5 times the max range requirement
rangeMax = 100;                             % Maximum range (m)
tm = 5*range2time(rangeMax,c)*512/500;  % Chirp duration (s) the sweep time should be at least 5 to 6 times the round trip
%range2time()这个函数是计算对于rangeMax处的目标，信号发出到接收需要的时间
% Determine the waveform bandwidth from the required range resolution
rangeRes = 1;                               % Desired range resolution (m)
bw = range2bw(rangeRes,c);                  % Corresponding bandwidth (Hz)

% Set the sampling rate to satisfy both the range and velocity requirements
% for the radar
sweepSlope = bw/tm;                           % FMCW sweep slope (Hz/s)
fbeatMax = range2beat(rangeMax,sweepSlope,c); % Maximum beat frequency (Hz)
%range2beat()这个函数能够根据距离算出差拍信号频率
vMax = 230*1000/3600;                    % Maximum Velocity of cars (m/s)
fdopMax = speed2dop(2*vMax,lambda);      % Maximum Doppler shift (Hz)

fifMax = fbeatMax+fdopMax;   % Maximum received IF (Hz)
fs = max(2*fifMax,bw);       % Sampling rate (Hz)
%% 目标参数
range1 = 40; %距离（米）此目标带干扰雷达
range2 = 70; %距离（米）
Velocity1 = 0*1000/3600;  %速度（米/秒）
Velocity2 = 0*1000/3600;  %速度（米/秒）
td1 = range2time(range1,c);
td2 = range2time(range2,c);
fdop1 = speed2dop(Velocity1,lambda);
fdop2 = speed2dop(Velocity2,lambda);

%% 恒虚警检测参数设置 guard,train,Pf
guard = 2;
train = 20;
Pf = 1e-2;

%% 本车雷达信号
L=fs*tm;
Nft=L;
rangeIdxToMeters = c * fs / (2 * sweepSlope * Nft);
t=(0:L-1)/fs;
% ms = cos(2*pi*sweepSlope*td1*t+2*pi*fc*td1-pi*sweepSlope*td1*td1)+...
%  0.1*cos(2*pi*sweepSlope*td2*t+2*pi*fc*td2-pi*sweepSlope*td2*td2);
A=100;
ms = A*cos(2*pi*(sweepSlope*td1-fdop1)*t+2*pi*(fc-bw/2+fdop1)*td1-pi*sweepSlope*td1*td1)+...
 A*0.5*cos(2*pi*(sweepSlope*td2-fdop2)*t+2*pi*(fc-bw/2+fdop2)*td2-pi*sweepSlope*td2*td2)+...
     A*0.03*randn(size(t));
% ms = cos(2*pi*(sweepSlope*td1-fdop1)*t+2*pi*(fc-bw/2+fdop1)*td1-pi*sweepSlope*td1*td1)+...
%  0.5*cos(2*pi*(sweepSlope*td2-fdop2)*t+2*pi*(fc-bw/2+fdop2)*td2-pi*sweepSlope*td2*td2);
Y=fft(ms, L);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
rnggrid = linspace(0, (L/2) * rangeIdxToMeters, L/2+1);%求距离维的坐标标度
[detected,th] = myCFAR(P1,guard,train,Pf);
detindx=find(detected);
sigA=sum(P1(detindx));
noiseA=sum(P1)-sigA;
SNR1=20*log10(P1(detindx)/noiseA);
formatSpec = "无干扰距离谱\nSNR1 is: %.3fdB SNR2 is: %.3fdB";
str = sprintf(formatSpec,SNR1(1),SNR1(2));

figure
plot(rnggrid,10*log10(P1),rnggrid,10*log10(th),'g--',...
    rnggrid(detindx),10*log10(P1(detindx)),'o')
title(str)
xlabel('距离 (米)')
ylabel('幅值(db)')
h=legend('距离谱','CFAR阈值','目标');
set(h, 'Box', 'off')
xlim([0,rangeMax]);
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))
%% 干扰信号建立
fcI = 77e9;
lambdaI = c/fcI;
bwI = bw;
tmI = tm/16;
sweepSlopeI = bwI/tmI;
tdI = range1/c;
fdopI = speed2dop(Velocity1,lambdaI);
% mISS = 10*cos(2*pi*((fc-fcI)-(bw/2-bwI/2)+(sweepSlopeI*tdI-fdopI))*t+pi*(sweepSlope-sweepSlopeI)*t.*t ...
%            +2*pi*(fcI-bwI/2+fdopI)*tdI-pi*sweepSlopeI*tdI*tdI);
AI=3200;
mIDS = AI*cos(2*pi*((fc-fcI)-(bw/2+bwI/2)-(sweepSlopeI*tdI+fdopI))*t+pi*(sweepSlope+sweepSlopeI)*t.*t ...
           +2*pi*(fcI+bwI/2+fdopI)*tdI+pi*sweepSlopeI*tdI*tdI);%randn(size(t)).*
L_int = round(L/(tm/tmI));
%% 混合信号分析
m = ms;
x=0:1:L_int-1;
y=gaussmf(x,[L_int/8 L_int/2])+3*A/AI;
% m(L_int*2+1:L_int*3)=10*mIDS(L_int*2+1:L_int*3).*y;%1
% m(L_int*4+3:L_int*5+2)=mIDS(L_int*4+3:L_int*5+2);% 
% m(L_int*6+3:L_int*7+2)=mIDS(L_int*6+3:L_int*7+2);
% m(L_int*8+2:L_int*9+1)=mIDS(L_int*8+2:L_int*9+1).*y;
% m(L_int*10+2:L_int*11+1)=mIDS(L_int*10+2:L_int*11+1);
% m(L_int*12+2:L_int*13+1)=mIDS(L_int*12+2:L_int*13+1).*y;
% m(L_int*4+1:L_int*5)=mIDS(L_int*4+1:L_int*5).*y;% 
m(L_int*6+1:L_int*7)=10*mIDS(L_int*6+1:L_int*7).*y;% 2
% m(L_int*8+1:L_int*9)=mIDS(L_int*8+1:L_int*9).*y;
% m(L_int*10+1:L_int*11)=mIDS(L_int*10+1:L_int*11).*y;
% m(L_int*12+1:L_int*13)=10*mIDS(L_int*12:L_int*13-1).*y;% 3
m_Original=m;

%% 时域分析
figure
plot(t*1e6,ms);
xlabel('time（ms）')
ylabel('Amplitude')
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))
figure
plot(t*1e6,m);

%% WEN法去噪
%找包络
up_Envelope = WEN(m,L_int/2);
% plot(t*1e6,up_Envelope,'r--');
xlabel('time（ms）')
ylabel('Amplitude')
% h=legend('Sinal','Envelope');
% set(h, 'Box', 'off')
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))
% hold off
%% 计算变点个数
% Y1 = diff(up_Envelope);
% Y2 = [0 Y1];

%% 归一化
WENsig=m./up_Envelope;
figure
plot(t*1e6,WENsig);
xlabel('时间（ms）')
ylabel('幅值')
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))
title('WEN归一化时域图');
Y=fft(WENsig, L);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
[detected,th] = myCFAR(P1,guard,train,Pf);
detindx=find(detected);
sigA=sum(P1(detindx));
noiseA=sum(P1)-sigA;
SNR2=20*log10(P1(detindx)/noiseA);
formatSpec = "WEN距离谱\nSNR1 is: %.3fdB SNR2 is: %.3fdB";
str = sprintf(formatSpec,SNR2);


figure
plot(rnggrid,10*log10(P1),rnggrid,10*log10(th),'g--',...
    rnggrid(detindx),10*log10(P1(detindx)),'o')
title(str)
xlabel('距离 (米)')
ylabel('幅值(db)')
h=legend('距离谱','CFAR阈值','目标');
set(h, 'Box', 'off')
xlim([0,rangeMax]);
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))

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
plot(t*1e6,AWENsig);
xlabel('时间（ms）')
ylabel('幅值')
title('AWEN归一化时域图');
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))
Y=fft(AWENsig, L);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% cfar = phased.CFARDetector('NumTrainingCells',40,'NumGuardCells',2);
% exp_pfa = 1e-3;
% cfar.ThresholdFactor = 'Auto';
% cfar.ProbabilityFalseAlarm = exp_pfa;
% cfar.ThresholdOutputPort = true;
% P1=P1.';
% CUTIdx=1:length(P1);
% [x_detected,th] = cfar(P1,1:length(P1));
[detected,th] = myCFAR(P1,guard,train,Pf);
detindx=find(detected);
sigA=sum(P1(detindx));
noiseA=sum(P1)-sigA;
SNR3=20*log10(P1(detindx)/noiseA);
formatSpec = "AWEN距离谱\nSNR1 is: %.3fdB SNR2 is: %.3fdB";
str = sprintf(formatSpec,SNR3);

figure
plot(rnggrid,10*log10(P1),rnggrid,10*log10(th),'g--',...
    rnggrid(detindx),10*log10(P1(detindx)),'o')
title(str)
xlabel('距离 (米)')
ylabel('幅值(db)')
h=legend('距离谱','CFAR阈值','目标');
set(h, 'Box', 'off')
xlim([0,rangeMax]);
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))
%% AR
%先测出干扰的起始与结束时间点
% 计算变点个数
beta = 0.4;%控制参数
Tc = (1+beta)*min(up_Envelope);
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
% Y2 = abs(Y2);
figure
plot(t*1e6,Y2)
xlabel('time（ms）')
ylabel('Amplitude')
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))
% beta=0.45;
% Tu=beta*min(Y2);
cptn=find(Y2);
% findchangepts(m,'Statistic','std','MaxNumChanges',length(cptn));
ipt = findchangepts(m,'Statistic','std','MaxNumChanges',length(cptn));
M = reshape(ipt,2,[]).';
maxl=M(:,1);
minl=M(:,2);
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
P = length(AWENsig);    % Number of samples in the extrapolated time series
% t = 1:P;
for k=1:length(minl)
for i=maxl(k):minl(k)
    m(i)=NaN;
end
end
lb = fillgaps(m,3000,200);

% lb=myfilgap(AWENsig,200,M,M2);

figure
subplot(2,1,1)
plot(t*1e6,m)
xlabel('时间（ms）')
ylabel('幅值')
title('Original');
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))
ipt=ipt';
xpos = transpose(double(ipt(1:end))*ones(1,2));
ypos = diag(ylim)*ones(2,numel(ipt));
line(t(xpos)*1e6,ypos,'Color','green');

subplot(2,1,2)
plot(t*1e6,lb)
line(t(xpos)*1e6,ypos,'Color','green');
xlabel('时间（ms）')
ylabel('幅值')
title('Reconstructed');
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))
% h=legend('Original','Reconstructed')
Y=fft(lb, L);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
[detected,th] = myCFAR(P1,guard,train,Pf);
detindx=find(detected);
sigA=sum(P1(detindx));
noiseA=sum(P1)-sigA;
SNR4=20*log10(P1(detindx)/noiseA);

figure
formatSpec = "AR距离谱\nSNR1 is: %.3fdB SNR2 is: %.3fdB";
str = sprintf(formatSpec,SNR4);
plot(rnggrid,10*log10(P1),rnggrid,10*log10(th),'g--',...
    rnggrid(detindx),10*log10(P1(detindx)),'o')
title(str)
xlabel('距离 (米)')
ylabel('幅值(dB)')
h=legend('距离谱','CFAR阈值','目标');
set(h, 'Box', 'off')
xlim([0,rangeMax]);
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))

%%
Y=fft(m_Original, L);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
sigA=sum(P1(detindx));
noiseA=sum(P1)-sigA;
SNR5=20*log10(P1(detindx)/noiseA);
formatSpec = "无处理距离谱\nSNR1 is: %.3fdB SNR2 is: %.3fdB";
str = sprintf(formatSpec,SNR5);
figure
plot(rnggrid,10*log10(P1),...
    rnggrid(detindx),10*log10(P1(detindx)),'o')
title(str)
xlabel('距离 (米)')
ylabel('幅值(dB)')
xlim([0,rangeMax]);
% set(gca,'FontSize',7)
set(gca,'LooseInset',get(gca,'TightInset'))
%% save
saveas(1,'1.jpg');
saveas(2,'2.jpg');
saveas(3,'3.jpg');
saveas(4,'4.jpg');
saveas(5,'5.jpg');
saveas(6,'6.jpg');
saveas(7,'7.jpg');
saveas(8,'8.jpg');
saveas(9,'9.jpg');
saveas(10,'10.jpg');
saveas(11,'11.jpg');