%% 生成模拟信号
fs = 1000;              % 采样率
t = 0:1/fs:1;           % 时间向量
f1 = 50;                % 信号频率1
f2 = 100;               % 信号频率2
f3 = 150;               % 信号频率3
x = 1*sin(2*pi*f1*t) + 2*sin(2*pi*f2*t) + 0.5*sin(2*pi*f3*t);   % 发送信号

plot(t,x);
xlabel('时间/s');
ylabel('功率/dB');
xlim([0, 0.5]);
title('发送信号');
%% 添加随机噪声
SNR_dB = 10;                  % 信噪比（dB）
% 信噪比越高，噪声含量越低，复合图像接近理想值
% 相反，信噪比越低，噪声含量越高，复合图像随机性更强
SNR = 10^(SNR_dB/10);         % 信噪比（线性值）
% 根据信噪比转换公式，SN（dB）=10lg（SN（倍数））
P_signal = rms(x)^2;         % 信号功率
% rms为有效值函数，功率等于有效值的平方
P_noise = P_signal / SNR;    % 噪声功率
% 根据信号功率和信噪比，算出噪声功率
noise = sqrt(P_noise) * randn(size(x));  % 生成随机噪声
rx = x + noise;               % 叠加噪声
figure
plot(t,rx);
xlabel('时间/s');
ylabel('功率/dB');
xlim([0, 0.5]);
title('接收信号');
%% 设计IIR滤波器
% 带通 IIR 滤波器，20 阶，下-3dB频率为80 Hz，上-3dB频率为120 Hz
% 滤波器类型为Butterworth，采样率为 1000 Hz
bpFilt = designfilt('bandpassiir','FilterOrder',20, ... % 滤波器阶数代表滤波器的传递函数中极点个数,过滤谐波次数
         'HalfPowerFrequency1',80, ... % 半功率点又名3dB点
         'HalfPowerFrequency2',120, ...
         'DesignMethod', 'butter', ... % 巴特沃斯滤波器具有平滑的单调频率响应，在通频带内最大程度地平坦。
         'SampleRate',fs);
fvtool(bpFilt) % 查看滤波器幅值响应
filt_rx = filter(bpFilt, rx); % 进行滤波
%% 绘制时域波形图
figure
subplot(3,1,1);
plot(t,x);
xlabel('时间/s');
ylabel('功率/dB');
title('发送信号');
subplot(3,1,2);
plot(t,rx);
xlabel('时间/s');
ylabel('功率/dB');
title('接收信号');
subplot(3,1,3);
plot(t,filt_rx);
xlabel('时间/s');
ylabel('功率/dB');
title('滤波信号');
% 由于有相位延迟，所以滤波信号一开始很平
%% 绘制频域波形图
figure
N = length(x); % 信号长度
f = (0:N-1)*(fs/N); % 频率向量

subplot(3,1,1);
X = fft(x)/N; % FFT变换并归一化
plot(f,abs(X)) % abs模值函数
xlabel('频率/Hz');
ylabel('功率/dB');
title('发送信号');

subplot(3,1,2);
RX = fft(rx)/N; % FFT变换并归一化
plot(f,abs(RX))
xlabel('频率/Hz');
ylabel('功率/dB');
title('接收信号');

subplot(3,1,3);
FILT_RX = fft(filt_rx)/N; % FFT变换并归一化
plot(f,abs(FILT_RX))
xlabel('频率/Hz');
ylabel('功率/dB');
title('滤波信号');
