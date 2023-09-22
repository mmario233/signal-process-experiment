%%  DSPʵ���-��������Ƶ�źŴ���
clear all; clc;
%% ��Ƶ���룬����ʱ���Σ�����
[x,fs]=audioread('yinpin.mp4');% ������Ƶ���ݣ�xΪ�����źţ�fsΪ������
x=x(:,1);                      % ȡ����x�ĵ�һ�и�ֵ��x��(��ȡ˫����Ƶ������һ·)  
len=length(x);
t=(0:len-1)/fs;                % ����ʱ�����t
figure                         % ����ԭʼ�����źŵ�ʱ����ͼ  
plot(t,x);                      
title('ԭʼ�����ź�ʱ����');   
xlabel('ʱ��/s');ylabel('����');grid on;axis tight;    
sound(x,fs);                   % ���������ź�
%%  �ڶ���������ź��м����˹������������ʱ���Σ�����
noise_mu=0;                    % ����������ֵ
noise_var=0.001;               % ������������
noise=randn(size(x)).*sqrt(noise_var)+noise_mu; % ������ԭʼ����������ͬ���������   
x_noise=x+noise;               % ��������ӵ�ԭʼ�����У��õ����������ź�  
figure                         % ���Ƽ���������������ź�ʱ��ͼ  
plot(t,x_noise);                    
title('���������ź�ʱ����');  
xlabel('ʱ��/s');ylabel('����');grid on;axis tight; 
sound(x_noise,fs);             % �������������ź�
%% ��ԭʼ�����źźͼ��������źŽ����׷���������Ƶ��ͼ�����бȽ�
f=fs*(-round(len)/2:len-round(len)/2-1)'/len;  
X=fft(x);                      % ��ԭʼ�����ź���FFT
X_noise=fft(x_noise);          % �Լ��������ź���FFT      
figure                         % ����Ƶ��ͼ                 
subplot(121);plot(f,fftshift(abs(X)));  
title('ԭʼ�����ź�Ƶ��'); 
xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight; 
subplot(122);plot(f,fftshift(abs(X_noise))); % �۲��ͼ��ȷ���˲������� 
title('���������ź�Ƶ��'); 
xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight; 
audiowrite('yinpin_noise.mp4',x_noise,fs); % �������������ź�(����)
%% ���IIR�����˲�����ȥ��
% �˲������
Filter_p=0.5*1e4;              % ����ͨ������Ƶ��/Hz
Filter_s=0.6*1e4;              % �����������Ƶ��/Hz
As=20;                         % �������/dB
Ap=1;                          % ͨ������/dB
wp=2*pi*Filter_p/fs;
ws=2*pi*Filter_s/fs; 
filter_p=2*fs*tan(wp/2);
filter_s=2*fs*tan(ws/2);
[n,wn]=buttord(wp,ws,Ap,As,'s');% ���ͨ�˲��������ͽ�ֹƵ��
[b,a]=butter(n,wn,'s');         % ��S��Ƶ����Ӧ���� 
[B,A]=bilinear(b,a,1);          % ˫���Ա任ʵ��S��Z��任 
[h,w]=freqz(B,A);               % ���ݲ�����Ƶ����Ӧ
figure
subplot(121)
plot(w*fs*0.5/pi,abs(h));       % �����˲���
title('IIR�˲���');
xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight;
subplot(122)
plot(w*fs*0.5/pi,20*log(abs(h)/abs(max(h))));       % �����˲���
title('IIR�˲���');
xlabel('Ƶ��/Hz');ylabel('����/dB');grid on;axis tight;
% �˲�
x_denoise=filter(B,A,x_noise);  % �����˲����Լ��������ź��˲� 
X_denoise=fft(x_denoise); 
figure  
subplot(2,2,1);plot(t,x_noise);
title('���������ź�ʱ����');
xlabel('ʱ��/s');ylabel('����');grid on;axis tight;
subplot(2,2,2);plot(t,x_denoise);
title('�˲��������ź�ʱ����');
xlabel('ʱ��/s');ylabel('����');grid on;axis tight;
subplot(2,2,3);plot(f,fftshift(abs(X_noise)));
title('���������ź�Ƶ��');
xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight; 
subplot(2,2,4);plot(f,fftshift(abs(X_denoise)));
title('�˲��������ź�Ƶ��');
xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight;
% ����������ȥ����Ƶ
sound(x_denoise,fs);                           % �����˲���������ź�         
audiowrite('yinpin_denoise.mp4',x_denoise,fs); % ����ȥ�������ź�(����)
%% ��ȡ
%%
for M = 2:10                             % ���ó�ȡ����
    % x_M=x(1:M:end);               % ��downsample��Ч
    x_M=downsample(x,M);            % ˼��������decimate�кβ�ͬ��
    X_M=fft(x_M);
    subplot(251);plot(f,fftshift(abs(X)));  
    title('ԭʼ�����ź�Ƶ��'); 
    xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight; 
    subplot(2,5,M);plot(f(1:M:end),fftshift(abs(X_M)));
    title({'��ȡ�����ź�Ƶ��',strcat('ÿ',num2str(M),"�����ȡһ����")});
    xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight; 
    sound(x_M,fs/M);                % ������ȡ�������ź�
    audiowrite(strcat('yinpin_downsampling',num2str(M),'.mp4'),x_M,fs); % ������ȡ�����ź�(����)
end
%% ��ֵ�ڲ�
%%
for L = 2:4  
% x_L=zeros(length(x)*L,1);x_L(1:L:end)=x;  % ��upsample��Ч
x_L=upsample(x,L);              % ˼��������interp�кβ�ͬ��
X_L=fft(x_L);
f_L=fs*(-round(length(x_L))/2:length(x_L)-round(length(x_L))/2-1)'/length(x_L);  
subplot(221);plot(f,fftshift(abs(X)));  
title('ԭʼ�����ź�Ƶ��'); 
xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight; 
subplot(2,2,L);plot(f_L,fftshift(abs(X_L)));
title(strcat('��ֵ�ڲ������ź�Ƶ��,ÿ����֮���',num2str(L-1),"����"));
xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight; 
sound(abs(x_L),fs*L);                % �����ڲ�������ź�
audiowrite(strcat('yinpin_upsampling',num2str(L),'.mp4'),x_L,fs); % �����ڲ������ź�(����)
end
%% �ײ�����
for chu = 5:-1:1
x_add=[zeros(len*chu,1);x];
X_add=fft(x_add);
x_add_1=ifft(X_add);
f_add=fs*(-round(length(x_add))/2:length(x_add)-round(length(x_add))/2-1)'/length(x_add);   
subplot(231);plot(f,fftshift(abs(X)));  
title('ԭʼ�����ź�Ƶ��'); 
xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight; 
subplot(2,3,chu+1);plot(f_add,fftshift(abs(X_add)));
title({'�ײ����������ź�Ƶ��',strcat('��ԭʼ�źų��ȵ�',num2str(chu),"��֮һ����")});
xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight; 
sound(x_add,fs);                % ��������������ź�
audiowrite('yinpin_addzero_top.mp4',x_add,fs); % �����ײ����������ź�(����)
end
%% β������
for chu = 5:-1:1
%x_add=[zeros(len*chu,1);x]; % �ײ�����
x_add=[x;zeros(len*chu,1)];% β������
% x_add=[zeros(round(len/2)*chu,1);x;zeros(round(len/2)*chu,1)]; ͬʱ����
X_add=fft(x_add);
x_add_1=ifft(X_add);
f_add=fs*(-round(length(x_add))/2:length(x_add)-round(length(x_add))/2-1)'/length(x_add);   
subplot(231);plot(f,fftshift(abs(X)));  
title('ԭʼ�����ź�Ƶ��'); 
xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight; 
subplot(2,3,chu+1);plot(f_add,fftshift(abs(X_add)));
title({'β�����������ź�Ƶ��',strcat('��ԭʼ�źų��ȵ�',num2str(chu),"������")});
xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight; 
sound(x_add,fs);                % ��������������ź�
audiowrite('yinpin_addzero_tail.mp4',x_add,fs); % ����β�����������ź�(����)
end
%% �ס�β������
for chu = 1:5
x_add=[zeros(round(len/2)*chu,1);x;zeros(round(len/2)*chu,1)];
X_add=fft(x_add);
x_add_1=ifft(X_add);
f_add=fs*(-round(length(x_add))/2:length(x_add)-round(length(x_add))/2-1)'/length(x_add);   
subplot(231);plot(f,fftshift(abs(X)));  
title('ԭʼ�����ź�Ƶ��'); 
xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight; 
subplot(2,3,chu+1);plot(f_add,fftshift(abs(X_add)));
title({'�ס�β�����������ź�Ƶ��',strcat('��ԭʼ�źų��ȵ�',num2str(chu),"������")});
xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight; 
sound(x_add,fs);                % ��������������ź�
audiowrite('yinpin_addzero_tt.mp4',x_add,fs); % �����ס�β�����������ź�(����)
end
%% �����ź��ڲ�
M=3;
x_M=downsample(x,M);            % ˼��������decimate�кβ�ͬ��
X_M=fft(x_M);
subplot(131);plot(f,fftshift(abs(X)));  
title('ԭʼ�����ź�Ƶ��'); 
xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight; 
subplot(132);plot(f(1:M:end),fftshift(abs(X_M)));
title({'��ȡ�����ź�Ƶ��',strcat('ÿ',num2str(M),"�����ȡһ����")});
xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight;
L=3;
x_L=interp(x,L);              % ˼��������interp�кβ�ͬ��
X_L=fft(x_L);
f_L=fs*(-round(length(x_L))/2:length(x_L)-round(length(x_L))/2-1)'/length(x_L);  
subplot(1,3,3);plot(f_L,fftshift(abs(X_L)));
title(strcat('��ֵ�ڲ������ź�Ƶ��,ÿ����֮���',num2str(L-1),"����"));
xlabel('Ƶ��/Hz');ylabel('����');grid on;axis tight; 
sound(abs(x_L),fs*L);                % �����ڲ�������ź�
audiowrite(strcat('yinpin_upsampling',num2str(L),'.mp4'),x_L,fs); % �����ڲ������ź�(����)
%% ˼����
% 1.	��ȡʱ����MATLAB�Դ�����downsample��Ϊdecimate�������ź�Ƶ���кβ�ͬ��
% 2.	��ֵ�ڲ�ʱ����MATLAB�Դ�����upsample��Ϊinterp�������ź�Ƶ���кβ�ͬ��
% 3.	�Գ�ȡ�ź�x_M������ֵ�ڲ壬�۲������ź�Ƶ�ף�
% 4.	�������Ϊʲô�ܹ���Сդ��ЧӦ��
%%
X1 = [664.2639804, 652.0547499, 628.3053846, 586.9887341];
X2 = [430.876117,323.2023371, 242.4356946, 150.1297219, 98.28940796, 64.74888664, ...
  31.99484137];
Y1 = [1533.546326, 1047.169811, 735.8178814, 522.1401658];
Y2 = [241.14,161.487965, 119.67, 74.74178404, 48.71900826, 31.41424273, ...
  14.04586241];
