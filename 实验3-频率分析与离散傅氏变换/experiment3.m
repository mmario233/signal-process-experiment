%% 生成并绘制2^n个子载波
clear
close all
n=6;
N=2^n; % N = 64
Fs = 200; % 采样频率
T = N/Fs; % 
num = 63; % 绘制子载波个数
y = repmat((1+1i)/sqrt(2),1,64); % 初始相位
y1(num, N) = 0;
% y1(k)=\sum y exp(j2\pi ik/N)
for k=0:num
    for n=0:N-1
        y1(k+1,n+1)=y(n+1)*exp(1i*2*pi*k*n/N);
    end
end
% y_f即子载波频带信息
for k=0:num
    y_f(k+1,:)=fft(y1(k+1,:));
end
% 补零获取更多信息
% 增加了对真实频谱采样的点数
% 并改变了采样点的位置
% 这将会显示出原信号频谱的更多的细节
% 不能提高信号的频率分辨率
a=20;
y11=[y1,zeros(num+1,a*64)];
f=(0:(a+1)*N-1)/T/(a+1)-Fs/2;
for k=0:num
    y_fft(k+1,:)=ifftshift(fft(y11(k+1,:)));
end
%绘图输入序列
figure
% 补零后输入信号的幅频特性
subplot(2,1,1)
hold on
for i = 1:num+1
    plot(f,abs(y_fft(i,:)));
end
hold off
ylabel('|X(f)|');
xlabel('f/Hz');
title('输入序列，频率间隔为200/64Hz');
% 补零后输入信号的相频特性
subplot(2,1,2)
hold on
for i = 1:num+1
    plot(f,angle(y_fft(i,:)));
end
hold off
ylabel('arg[X(f)]/rad');
xlabel('f/Hz');
%% 绘制IFFT
for k=0:num
    my2(k+1,:)=myifft(y_f(k+1,:));
end
figure
subplot(2,1,1)
hold on
for i = 1:5
plot(0:1/Fs:(N-1)/Fs,real(my2(i,:)));
end
hold off
ylabel('x(t)');
xlabel('t');
title('自拟IFFT函数结果(实部)');
subplot(2,1,2)
hold on
for i = 1:5
plot(0:1/Fs:(N-1)/Fs,imag(my2(i,:)));
end
hold off
ylabel('x(t)');
xlabel('t');
title('自拟IFFT函数结果(虚部)');
%% 绘制FFT函数
f=(0:N-1)/T-Fs/2;
for k=0:num
    my3(k+1,:)=myfft(my2(k+1,:));
end
figure
subplot(2,1,1)
hold on
for i = 1:num+1
plot(f,abs(my3(i,:)));
end
hold off
ylabel('|X(f)|');
xlabel('f/Hz');
title('自拟FFT函数结果');
subplot(2,1,2)
hold on
for i = 1:num+1
plot(f,angle(my3(i,:)));
end
hold off
ylabel('arg[X(f)]/rad');
xlabel('f/Hz');
%% 比较MATLAB自带FFT与自拟FFT、DFT的性能
profile clear
profile on
for k=0:num
    z1(k+1,:)=fft(my2(k+1,:));
end
for k=0:num
    z2(k+1,:)=myfft(my2(k+1,:));
end
for k=0:num
    z3(k+1,:)=mydft(my2(k+1,:));
end
profile report
%% 结果比较
z1-z2
%% 1024
clear
profile clear
profile on
n=10;
N=2^n; % N = 1024
Fs = 3000; % 采样频率
T = N/Fs; % 
num = 1023; % 绘制子载波个数
t=0:N-1;
y = repmat(1i,1,1024); % 初始相位
y1(num, N) = 0;
% y1(k)=\sum y exp(j2\pi ik/N)
for k=0:num
    for n=0:N-1
        y1(k+1,n+1)=y(n+1)*exp(1i*2*pi*k*n/N);
    end
end
for k=0:num
    z1(k+1,:)=fft(y1(k+1,:));
end
for k=0:num
    z2(k+1,:)=myfft(y1(k+1,:));
end
for k=0:num
    z3(k+1,:)=mydft(y1(k+1,:));
end
profile viewer