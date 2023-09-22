% 基本原理
% 经参数可调滤波器，与期望信号比对产生误差信号，
% 通过自适应算法使误差信号的均方最小
% LMS自适应算法直接利用瞬态均方误差对滤波器系数求梯度
% 特点
% 没有关于待提取信息的先验统计知识
% 直接利用观测数据依据某种判据在观测过程中不断递归更新
% 最优化
% 方法
% 梯度下降法
% 优点
% 算法简单，易于实现，算法复杂度低，能够抑制旁瓣效应
% 缺点
% 收敛速率较慢，因为LMS滤波器系数更新是逐点的
% 跟踪性能较差，并且随着滤波器阶数(步长参数)升高，系统的稳定性下降
% 要求不同时刻的输入向量线性无关
clear;clc;echo off;close all;
N=10000;                 %指定信号序列长度
info=random_binary(N);   %产生双极性不归零基带信号序列（等概率随机正负1）
SNR_in_dB=8:1:18;        %AWGN信道信噪比 AWGN信道（添加了白噪声的信道）
for j=1:length(SNR_in_dB)
    [y,len]=channel(info,SNR_in_dB(j));  %通过既有码间干扰又有高斯白噪声信道
    numoferr=0;                          %初始误码统计数
    for i=len+1:N+len                   %从第len个码元开始为真实信号码元
        if (y(i)<0)                     %判决译码
            decis=-1;
        else
            decis=1;
        end
        if(decis~=info(i-len))          %判断是否误码，统计误码码元个数
            numoferr=numoferr+1;
        end
    end
    Pe(j)=numoferr/N;                    %未经均衡器均衡，得到的误码率
end
figure
semilogy(SNR_in_dB,Pe,'red*-');          %未经均衡器，误码率结果图
    hold on;                             %semilogy表示y坐标轴是对数坐标系
delta_1=0.11;     %指定自适应均衡器的步长
delta_2=0.2;     %指定自适应均衡器的步长
delta_3=0.09;     %指定自适应均衡器的步长
 
for j=1:length(SNR_in_dB)
    y=channel(info,SNR_in_dB(j));        %通过信道
    z=lms_equalizer(y,info,0.001);     %通过自适应均衡器，并设置步长为0.11
    numoferr=0;
    for i=1:N
        if (z(i)<0)
            decis=-1;
        else
            decis=1;
        end
        if (decis~=info(i))
            numoferr=numoferr+1;
        end
    end
    Pee(j)=numoferr/N;                   %经自适应均衡器均衡后，得到的误码率
end
% SNR = EbN0 + 10log10(nBits*coderate) - 10log10(0.5or1 * upfactor);
semilogy(SNR_in_dB,Pee,'black*-');       %自适应均衡器均衡之后，误码率结果图
    hold on;
xlabel('SNR in dB');
ylabel('Pe');
title('ISI信道自适应均衡系统仿真');
legend('未经均衡器均衡','经自适应均衡器均衡，步长detla=0.11');
% 随着步长的减小，误码率的起始点就比较小，随着信道信噪比的增加，误码率快速降低
% 0.2步长太大，反而提升误码率，且基本不随信噪比变化而变化
eyediagram(y(500:1000),10);             %均衡前眼图
eyediagram(z(500:1000),10);             %均衡后眼图，步长0.11
% 眼图的大小代表了信号的质量
% 将多个周期的波形叠加在同一个周期内显示
% 眼图张开的宽度决定了接收波形可以不受串扰影响而抽样再生的时间间隔
% 最佳抽样时刻在眼睛张开最大的时刻
% 斜边斜率表示系统对误差的灵敏度，斜率越大，越灵敏
% 门限电平为0
% 眼图左右角阴影部分的水平宽度表示信号零点的变化范围
% 抽样时刻 阴影区垂直宽度表示最大信号失真量
% 抽样时刻 上下阴影区间隔的一半是最小噪声容限