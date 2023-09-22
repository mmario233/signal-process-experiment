%%  DSP实验课-语音及音频信号处理
clear all; clc;
%% 音频导入，绘制时域波形，发声
[x,fs]=audioread('yinpin.mp4');% 导入音频数据，x为语音信号，fs为采样率
x=x(:,1);                      % 取矩阵x的第一列赋值到x中(提取双轨音频的其中一路)  
len=length(x);
t=(0:len-1)/fs;                % 生成时间变量t
figure                         % 绘制原始语音信号的时域波形图  
plot(t,x);                      
title('原始语音信号时域波形');   
xlabel('时间/s');ylabel('幅度');grid on;axis tight;    
sound(x,fs);                   % 发声语音信号
%%  在读入的语音信号中加入高斯白噪声，绘制时域波形，发声
noise_mu=0;                    % 设置噪声均值
noise_var=0.001;               % 设置噪声方差
noise=randn(size(x)).*sqrt(noise_var)+noise_mu; % 产生与原始语音长度相同的随机噪声   
x_noise=x+noise;               % 把噪声添加到原始语音中，得到加噪语音信号  
figure                         % 绘制加随机噪声后语音信号时域图  
plot(t,x_noise);                    
title('加噪语音信号时域波形');  
xlabel('时间/s');ylabel('幅度');grid on;axis tight; 
sound(x_noise,fs);             % 发声加噪语音信号
%% 对原始语音信号和加噪语音信号进行谱分析，绘制频谱图并进行比较
f=fs*(-round(len)/2:len-round(len)/2-1)'/len;  
X=fft(x);                      % 对原始语音信号做FFT
X_noise=fft(x_noise);          % 对加噪语音信号做FFT      
figure                         % 绘制频谱图                 
subplot(121);plot(f,fftshift(abs(X)));  
title('原始语音信号频谱'); 
xlabel('频率/Hz');ylabel('幅度');grid on;axis tight; 
subplot(122);plot(f,fftshift(abs(X_noise))); % 观察此图，确定滤波器参数 
title('加噪语音信号频谱'); 
xlabel('频率/Hz');ylabel('幅度');grid on;axis tight; 
audiowrite('yinpin_noise.mp4',x_noise,fs); % 导出加噪语音信号(单轨)
%% 设计IIR数字滤波器，去噪
% 滤波器设计
Filter_p=0.5*1e4;              % 设置通带截至频率/Hz
Filter_s=0.6*1e4;              % 设置阻带截至频率/Hz
As=20;                         % 阻带波纹/dB
Ap=1;                          % 通带波纹/dB
wp=2*pi*Filter_p/fs;
ws=2*pi*Filter_s/fs; 
filter_p=2*fs*tan(wp/2);
filter_s=2*fs*tan(ws/2);
[n,wn]=buttord(wp,ws,Ap,As,'s');% 求低通滤波器阶数和截止频率
[b,a]=butter(n,wn,'s');         % 求S域频率响应参数 
[B,A]=bilinear(b,a,1);          % 双线性变换实现S域到Z域变换 
[h,w]=freqz(B,A);               % 根据参数求频率响应
figure
subplot(121)
plot(w*fs*0.5/pi,abs(h));       % 绘制滤波器
title('IIR滤波器');
xlabel('频率/Hz');ylabel('幅度');grid on;axis tight;
subplot(122)
plot(w*fs*0.5/pi,20*log(abs(h)/abs(max(h))));       % 绘制滤波器
title('IIR滤波器');
xlabel('频率/Hz');ylabel('幅度/dB');grid on;axis tight;
% 滤波
x_denoise=filter(B,A,x_noise);  % 利用滤波器对加噪语音信号滤波 
X_denoise=fft(x_denoise); 
figure  
subplot(2,2,1);plot(t,x_noise);
title('加噪语音信号时域波形');
xlabel('时间/s');ylabel('幅度');grid on;axis tight;
subplot(2,2,2);plot(t,x_denoise);
title('滤波后语音信号时域波形');
xlabel('时间/s');ylabel('幅度');grid on;axis tight;
subplot(2,2,3);plot(f,fftshift(abs(X_noise)));
title('加噪语音信号频谱');
xlabel('频率/Hz');ylabel('幅度');grid on;axis tight; 
subplot(2,2,4);plot(f,fftshift(abs(X_denoise)));
title('滤波后语音信号频谱');
xlabel('频率/Hz');ylabel('幅度');grid on;axis tight;
% 发声，导出去噪音频
sound(x_denoise,fs);                           % 发声滤波后的语音信号         
audiowrite('yinpin_denoise.mp4',x_denoise,fs); % 导出去噪语音信号(单轨)
%% 抽取
%%
for M = 2:10                             % 设置抽取倍数
    % x_M=x(1:M:end);               % 与downsample等效
    x_M=downsample(x,M);            % 思考：换成decimate有何不同？
    X_M=fft(x_M);
    subplot(251);plot(f,fftshift(abs(X)));  
    title('原始语音信号频谱'); 
    xlabel('频率/Hz');ylabel('幅度');grid on;axis tight; 
    subplot(2,5,M);plot(f(1:M:end),fftshift(abs(X_M)));
    title({'抽取语音信号频谱',strcat('每',num2str(M),"个点抽取一个点")});
    xlabel('频率/Hz');ylabel('幅度');grid on;axis tight; 
    sound(x_M,fs/M);                % 发声抽取的语音信号
    audiowrite(strcat('yinpin_downsampling',num2str(M),'.mp4'),x_M,fs); % 导出抽取语音信号(单轨)
end
%% 零值内插
%%
for L = 2:4  
% x_L=zeros(length(x)*L,1);x_L(1:L:end)=x;  % 与upsample等效
x_L=upsample(x,L);              % 思考：换成interp有何不同？
X_L=fft(x_L);
f_L=fs*(-round(length(x_L))/2:length(x_L)-round(length(x_L))/2-1)'/length(x_L);  
subplot(221);plot(f,fftshift(abs(X)));  
title('原始语音信号频谱'); 
xlabel('频率/Hz');ylabel('幅度');grid on;axis tight; 
subplot(2,2,L);plot(f_L,fftshift(abs(X_L)));
title(strcat('零值内插语音信号频谱,每两点之间插',num2str(L-1),"个点"));
xlabel('频率/Hz');ylabel('幅度');grid on;axis tight; 
sound(abs(x_L),fs*L);                % 发声内插的语音信号
audiowrite(strcat('yinpin_upsampling',num2str(L),'.mp4'),x_L,fs); % 导出内插语音信号(单轨)
end
%% 首部补零
for chu = 5:-1:1
x_add=[zeros(len*chu,1);x];
X_add=fft(x_add);
x_add_1=ifft(X_add);
f_add=fs*(-round(length(x_add))/2:length(x_add)-round(length(x_add))/2-1)'/length(x_add);   
subplot(231);plot(f,fftshift(abs(X)));  
title('原始语音信号频谱'); 
xlabel('频率/Hz');ylabel('幅度');grid on;axis tight; 
subplot(2,3,chu+1);plot(f_add,fftshift(abs(X_add)));
title({'首部补零语音信号频谱',strcat('补原始信号长度的',num2str(chu),"分之一个零")});
xlabel('频率/Hz');ylabel('幅度');grid on;axis tight; 
sound(x_add,fs);                % 发声补零的语音信号
audiowrite('yinpin_addzero_top.mp4',x_add,fs); % 导出首部补零语音信号(单轨)
end
%% 尾部补零
for chu = 5:-1:1
%x_add=[zeros(len*chu,1);x]; % 首部补零
x_add=[x;zeros(len*chu,1)];% 尾部补零
% x_add=[zeros(round(len/2)*chu,1);x;zeros(round(len/2)*chu,1)]; 同时补零
X_add=fft(x_add);
x_add_1=ifft(X_add);
f_add=fs*(-round(length(x_add))/2:length(x_add)-round(length(x_add))/2-1)'/length(x_add);   
subplot(231);plot(f,fftshift(abs(X)));  
title('原始语音信号频谱'); 
xlabel('频率/Hz');ylabel('幅度');grid on;axis tight; 
subplot(2,3,chu+1);plot(f_add,fftshift(abs(X_add)));
title({'尾部补零语音信号频谱',strcat('补原始信号长度的',num2str(chu),"倍个零")});
xlabel('频率/Hz');ylabel('幅度');grid on;axis tight; 
sound(x_add,fs);                % 发声补零的语音信号
audiowrite('yinpin_addzero_tail.mp4',x_add,fs); % 导出尾部补零语音信号(单轨)
end
%% 首、尾部补零
for chu = 1:5
x_add=[zeros(round(len/2)*chu,1);x;zeros(round(len/2)*chu,1)];
X_add=fft(x_add);
x_add_1=ifft(X_add);
f_add=fs*(-round(length(x_add))/2:length(x_add)-round(length(x_add))/2-1)'/length(x_add);   
subplot(231);plot(f,fftshift(abs(X)));  
title('原始语音信号频谱'); 
xlabel('频率/Hz');ylabel('幅度');grid on;axis tight; 
subplot(2,3,chu+1);plot(f_add,fftshift(abs(X_add)));
title({'首、尾部补零语音信号频谱',strcat('补原始信号长度的',num2str(chu),"倍个零")});
xlabel('频率/Hz');ylabel('幅度');grid on;axis tight; 
sound(x_add,fs);                % 发声补零的语音信号
audiowrite('yinpin_addzero_tt.mp4',x_add,fs); % 导出首、尾部补零语音信号(单轨)
end
%% 抽样信号内插
M=3;
x_M=downsample(x,M);            % 思考：换成decimate有何不同？
X_M=fft(x_M);
subplot(131);plot(f,fftshift(abs(X)));  
title('原始语音信号频谱'); 
xlabel('频率/Hz');ylabel('幅度');grid on;axis tight; 
subplot(132);plot(f(1:M:end),fftshift(abs(X_M)));
title({'抽取语音信号频谱',strcat('每',num2str(M),"个点抽取一个点")});
xlabel('频率/Hz');ylabel('幅度');grid on;axis tight;
L=3;
x_L=interp(x,L);              % 思考：换成interp有何不同？
X_L=fft(x_L);
f_L=fs*(-round(length(x_L))/2:length(x_L)-round(length(x_L))/2-1)'/length(x_L);  
subplot(1,3,3);plot(f_L,fftshift(abs(X_L)));
title(strcat('零值内插语音信号频谱,每两点之间插',num2str(L-1),"个点"));
xlabel('频率/Hz');ylabel('幅度');grid on;axis tight; 
sound(abs(x_L),fs*L);                % 发声内插的语音信号
audiowrite(strcat('yinpin_upsampling',num2str(L),'.mp4'),x_L,fs); % 导出内插语音信号(单轨)
%% 思考题
% 1.	抽取时，将MATLAB自带函数downsample换为decimate，所得信号频谱有何不同？
% 2.	零值内插时，将MATLAB自带函数upsample换为interp，所得信号频谱有何不同？
% 3.	对抽取信号x_M进行零值内插，观察所得信号频谱？
% 4.	补零操作为什么能够减小栅栏效应？
%%
X1 = [664.2639804, 652.0547499, 628.3053846, 586.9887341];
X2 = [430.876117,323.2023371, 242.4356946, 150.1297219, 98.28940796, 64.74888664, ...
  31.99484137];
Y1 = [1533.546326, 1047.169811, 735.8178814, 522.1401658];
Y2 = [241.14,161.487965, 119.67, 74.74178404, 48.71900826, 31.41424273, ...
  14.04586241];
