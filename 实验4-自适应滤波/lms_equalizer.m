function [ z ] = lms_equalizer( y,info,delta )
%LMS算法自适应滤波器实现
estimated_c=[0 0 0 0 0 1 0 0 0 0 0]; %初始抽头系数(长度应该是和channel中信道阶数相同=11)
K=5;                                 %K=（length（estimated_c）-1）/2
for k=1:size(y,2)-2*K              %channel中返回参数len的长度也是5，或许K的选择便是基于len，需要K=len
     y_k=y(k:k+2*K);                 %获取码元，一次11个
     z_k=estimated_c*y_k';           %各抽头系数与码元相乘后求和
     e_k=info(k)-z_k;                %误差估计
     estimated_c=estimated_c+delta * e_k * y_k;%计算校正抽头系数
     z(k)=z_k;                       %均衡后输出的码元序列
end
%误差e=d-y，（y指经过信道后的输出信号）这里期望信号d使用的是输入信号info
%比较误码率都是只比较info的长度N
%size(y,2)返回y的列数=N+L-1（L=2len+1）,则size(y,2)-2*K=N+L-(2K+1)=N+2(len-K)
%这里将会有个弊端：如果len>K，则k=1：N+2(len-K)会超过info的长度N.
end