function [ y,len ] = channel( x,snr_in_dB )
%模拟既有码间干扰又有高斯白噪声的信道
SNR=exp(snr_in_dB*log(10)/10);   %信噪比真值转换
sigma=1/sqrt(2*SNR);             %高斯白噪声的标准差
%指定信道的ISI参数，可以看出此信道质量还是比较差的
% 时延大于持续时间，来自于不同路径的不同信号之间会互相干扰
actual_isi=[0.05 -0.063 0.088 -0.126 -0.25 0.9047 0.25 0 0.126 0.038 0.088];
%actual_isi(1) 到 actual_isi(6) 存储了码间串扰对应的 ISI 值
% ，而 actual_isi(7) 和 actual_isi(9) 分别表示 DC 偏移和基带增益失真对应的 ISI 值，
% actual_isi(8) 即为直流失真。另外，actual_isi(10) 和 actual_isi(11) 表示远端串扰和近端串扰对应的 ISI 值。
%第五个元素（-0.25）和第六个元素（0.9047）的绝对值较大，同时符号相反，这表示在该信道上存在比较严重的波形畸变，导致信号品质不佳
len_actual_isi=(length(actual_isi)-1)/2;
len=len_actual_isi;
y=conv(actual_isi,x);             %信号通过信道，相当于信号序列与信道模型序列作卷积
%需要指出，此时码元序列长度变为N+L=N+2len+1-1，译码时我们从第len个码元开始到N+len个结束
for i=1:2:size(y,2)
    [noise(i),noise(i+1)]=gngauss(sigma); %产生噪声
end
y=y+noise;                                %叠加噪声
%也可直接用y = awgn(y,SNR)
% 向信号"y"添加高斯白噪声，信噪比大小为“SNR”，单位是dB；信号“y”的功率假定为 0 dBW；
end

%所以生成x是10000的列向量，用ISI和信道中的序列做了卷积，译码选取的是第5个码元到第N+5 ； 