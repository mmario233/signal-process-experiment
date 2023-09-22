clear
close
n=6;
N=2^n;
Fs = 200;
T = N/Fs;
num = 63;
t=0:N-1;
y = repmat(1i,1,64); % 初始相位
y1(num, N) = 0;
% y1(k)=\sum y exp(j2\pi ik/N)
for k=0:num
for n=0:N-1
    y1(k+1,n+1)=y(n+1)*exp(1i*2*pi*k*n/N);
end
end
for k=0:num
    y_f(k+1,:)=fft(y1(k+1,:));
end
a=20;
y11=[y1,zeros(num+1,a*64)];
f=(0:(a+1)*N-1)/T/(a+1)-Fs/2;
for k=0:num
    y_fft(k+1,:)=fftshift(fft(y11(k+1,:)));
end

%绘图输入序列
figure
subplot(2,1,1)
hold on
for i = 1:num+1
    plot(f,abs(y_fft(i,:)));
end
hold off
ylabel('|X(f)|');
xlabel('f/Hz');
title('输入序列，频率间隔为200/64Hz');
subplot(2,1,2)
hold on
for i = 1:num+1
    plot(f,angle(y_fft(i,:)));
end
hold off
ylabel('arg[X(f)]/rad');
xlabel('f/Hz');
%%
%绘图IFFT
for k=0:num
    y2(k+1,:)=ifft(y_f(k+1,:));
end
my2 = zeros(size(y2));
for k=0:num
    my2(k+1,:)=myifft(y_f(k+1,:));
end
figure
subplot(2,1,1)
hold on
for i = 1:10
plot(0:1/Fs:(N-1)/Fs,y2(i,:));
end
hold off
ylabel('x(t)');
xlabel('t');
title('MATLAB自带函数结果');
subplot(2,1,2)
hold on
for i = 1:num+1
plot(0:1/Fs:(N-1)/Fs,my2(i,:));
end
hold off
ylabel('x(t)');
xlabel('t');
title('自拟IFFT函数结果');
%%
%绘图FFT
f=(0:N-1)/T-Fs/2;
for k=0:num
    y3(k+1,:)=fft(y2(k+1,:));
end
for k=0:num
    my3(k+1,:)=myfft(y2(k+1,:));
end
figure
subplot(2,2,1)
hold on
for i = 1:num+1
plot(f,abs(y3(i,:)));
end
hold off
ylabel('|X(f)|');
xlabel('f/Hz');
title('MATLAB自带函数结果');
subplot(2,2,2)
hold on
for i = 1:num+1
plot(f,angle(y3(i,:)));
end
hold off
ylabel('arg[X(f)]/rad');
xlabel('f/Hz');
subplot(2,2,3)
hold on
for i = 1:num+1
plot(f,abs(my3(i,:)));
end
hold off
ylabel('|X(f)|');
xlabel('f/Hz');
title('自拟FFT函数结果');
subplot(2,2,4)
hold on
for i = 1:num+1
plot(f,angle(my3(i,:)));
end
hold off
ylabel('arg[X(f)]/rad');
xlabel('f/Hz');
%% MATLAB自带IFFT函数与自拟函数误差
y2-my2
%% MATLAB自带FFT/DFT函数与自拟函数误差
y3-my3